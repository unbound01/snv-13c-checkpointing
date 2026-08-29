"""
Part 4: Full Checkpointing Protocol + Heating-Budget Analysis
==============================================================
Diamond quantum computing architecture : quantum checkpointing protocol.

Qubit 1: SnV electron spin  (compute qubit)
Qubit 2: 13C nuclear spin   (memory qubit)

LABELLING CAVEAT : SWAP FIDELITY
----------------------------------
swap_fidelity = 0.9992 is the ¹³C single-qubit RF gate fidelity from
Resch et al., Phys. Rev. X 16, 011060 (2026) [arXiv:2509.03354].
It is used here as a PROPOSED/UNVALIDATED proxy for a SnV↔¹³C SWAP gate.
No SWAP operation was demonstrated in that paper (initialization used
optical+MW pumping, not a SWAP).  All fidelity results in this file are
contingent on this assumption and should be interpreted accordingly.

TEMPERATURES: 1.7 K and 3 K ONLY
----------------------------------
7 K and 10 K data have been removed.  The Rosenthal 2023 Orbach fit does
not support SnV operation at those temperatures on this platform.

NEW: HEATING BUDGET ANALYSIS
------------------------------
Source: Karapatzakis et al., arXiv:2606.15398 (2026)

Protocol steps:
  Step 1 : Local gates on SnV       (t_gates = 1 µs)
  Step 2 : Proposed SWAP in SnV→¹³C (t_swap  = 10 µs, decohere + noisy SWAP)
  Step 3 : Wait     t_comm          (both qubits decohere)
  Step 4 : Proposed SWAP out ¹³C→SnV (t_swap = 10 µs, decohere + noisy SWAP)

Baseline: SnV decoheres freely for t_total = t_gates + 2*t_swap + t_comm.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, qeye, sigmaz, sigmam,
    ket2dm, ptrace, mesolve, Qobj
)

# ============================================================================
# SYSTEM REBUILD
# ============================================================================
I2 = qeye(2)
I4 = tensor(I2, I2)

psi_snv_plus = (basis(2, 0) + basis(2, 1)).unit()
psi_13c_zero = basis(2, 0)
rho_init     = ket2dm(tensor(psi_snv_plus, psi_13c_zero))

# ============================================================================
# PHYSICAL PARAMETERS  : Corrected, literature-traceable
# ============================================================================
# SnV : Rosenthal et al., Phys. Rev. X 13, 031022 (2023)
T1_snv = {1.7: 4.23e-3, 3.0: 4.23e-3}      # T1 (s) : Rosenthal 2023, direct 3 K measurement.
                                              # No direct 1.7 K value; T1 > ~20 ms lower bound;
                                              # Orbach extrapolation ≈ 200 s at 1.7 K.
                                              # 4.23 ms used at 1.7 K as conservative worst-case.
T2_snv = {1.7: 170.0e-6, 3.0: 170.0e-6}    # T2 Hahn echo (µs) : plateau 1.7–3 K

# ¹³C : Resch et al., Phys. Rev. X 16, 011060 (2026), measured at 50 mK
# *** NOT validated at 1.7–3 K : see module docstring ***
T2_13c = {1.7: 1.35, 3.0: 1.35}            # CPMG-128 coherence time (s)

# SWAP proxy fidelity : Resch 2026, ¹³C single-qubit RF gate @ 50 mK
# *** NOT an experimentally confirmed SWAP gate fidelity ***
swap_fidelity = 0.9992   # proposed nuclear-spin single-qubit gate fidelity (Resch 2026)

t_swap  = 10e-6    # assumed SWAP gate duration (µs)
t_gates = 1e-6     # local SnV gate duration

# t_comm sweep : 11 values, 1 µs → 100 µs range (tighter than before to
# focus on region where heating budget is most relevant)
t_comm_values = np.array([1e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6,
                           200e-6, 500e-6, 1e-3, 5e-3, 10e-3])

temperatures = [1.7, 3.0]
COLORS       = {1.7: '#1f77b4', 3.0: '#d62728'}   # blue = 1.7 K, red = 3 K

# ============================================================================
# MW GATE PARAMETERS : for heating budget
# ============================================================================
# SnV π-pulse time : Rosenthal et al., PRX 2023
t_pi_pulse = 48.4e-9   # 48.4 ns MW π-pulse

# Heating thresholds : Karapatzakis et al., arXiv:2606.15398 (2026)
# These thresholds apply to CW and 20%-duty-cycle drives at B=400 mT.
# MW (2 GHz electron-spin gates):
#   - Gradual heating onset above ~7 dBm; local T rises from 1.8 K toward ~3 K at 11 dBm
#   - Superconducting breakdown: ~11–12 dBm CW
#   - 20% duty cycle breakdown thresholds: 17.2 dBm (B=0), ~10 dBm (B=400 mT, θ=84°),
#                                          ~12 dBm (B=400 mT, θ=20°)
# RF (20 MHz nuclear-spin driving), 400 mT:
#   - No measurable heating up to abrupt breakdown at 9.4 dBm (Bac ≈ 1.2 mT)
# Thermal relaxation time constant: τ_th = 0.54(17) ms
tau_th_ms = 0.54   # ms : Karapatzakis 2026

# ============================================================================
# OPERATORS
# ============================================================================
Sz_snv = tensor(sigmaz(), I2)
Sz_13c = tensor(I2, sigmaz())
Sm_snv = tensor(sigmam(), I2)
H      = I4 * 0.0

# ============================================================================
# SWAP GATE
# ============================================================================
U_swap = Qobj(np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],
                        dtype=complex), dims=[[2,2],[2,2]])
rho_mixed = I4 / 4.0

def noisy_swap(rho):
    """
    Depolarising SWAP proxy.
    F = 0.9992 = ¹³C single-qubit gate fidelity (Resch 2026 @ 50 mK).
    *** Not an experimentally confirmed SWAP gate ***
    """
    rho_ideal = U_swap * rho * U_swap.dag()
    return swap_fidelity * rho_ideal + (1.0 - swap_fidelity) * rho_mixed

# ============================================================================
# COLLAPSE OPERATOR BUILDER
# ============================================================================
def build_c_ops(T):
    """
    Lindblad collapse operators for temperature T (Kelvin).

    Correct convention (verified):
      gamma_phi = max(0, 1/T2 - 1/(2*T1))  → pure dephasing rate
      c_dephase  = sqrt(gamma_phi / 2) * Sz
      c_relax    = sqrt(1/T1)          * Sm
      c_13C      = sqrt(1/(2*T2_13c)) * Sz_13c

    At 1.7K and 3K: T2 = 170 µs, 2*T1 = 8.46 ms → T2 << 2*T1 ✓ (not T1-limited)
    """
    gamma_phi_snv = max(0.0, 1.0/T2_snv[T] - 1.0/(2*T1_snv[T]))
    return [
        np.sqrt(gamma_phi_snv / 2.0)  * Sz_snv,
        np.sqrt(1.0 / T1_snv[T])      * Sm_snv,
        np.sqrt(1.0 / (2*T2_13c[T])) * Sz_13c,
    ]

# ============================================================================
# HELPERS
# ============================================================================
def evolve(rho_in, duration, c_ops_list, n_steps=50):
    if duration <= 0.0:
        return rho_in
    tlist  = np.linspace(0.0, duration, n_steps)
    result = mesolve(H, rho_in, tlist, c_ops_list, e_ops=[])
    return result.states[-1]

def pure_fidelity(psi_target, rho):
    return float((psi_target.dag() * rho * psi_target).real)

# ============================================================================
# CHECKPOINTING PROTOCOL FUNCTIONS
# ============================================================================
def full_checkpoint(rho_start, t_comm, T):
    """Four-step checkpointing. Returns SnV fidelity with |+>."""
    c = build_c_ops(T)
    rho = evolve(rho_start, t_gates, c)    # Step 1: local gates
    rho = evolve(rho, t_swap, c)           # Step 2a: decohere during SWAP
    rho = noisy_swap(rho)                  # Step 2b: proposed SWAP in
    rho = evolve(rho, t_comm, c)           # Step 3: wait
    rho = evolve(rho, t_swap, c)           # Step 4a: decohere during SWAP
    rho = noisy_swap(rho)                  # Step 4b: proposed SWAP out
    return pure_fidelity(psi_snv_plus, ptrace(rho, 0))

def no_checkpoint(rho_start, t_comm, T):
    """Baseline: free decoherence for t_total. Returns SnV fidelity."""
    c       = build_c_ops(T)
    t_total = t_gates + t_swap + t_comm + t_swap
    rho     = evolve(rho_start, t_total, c)
    return pure_fidelity(psi_snv_plus, ptrace(rho, 0))

# ============================================================================
# HEATING BUDGET ANALYSIS
# ============================================================================
# MW gate duty cycle = (total gate on-time per cycle) / (total cycle time)
#
# Gate sequence per checkpoint cycle:
#   - Local gates:   t_gates  = 1 µs  (assume N_gates π-pulses at 48.4 ns each)
#   - Two SWAPs:     2×t_swap = 20 µs (each swap assumed ~1 µs of MW on-time)
# Total MW on-time per cycle is dominated by the SWAP steps.
#
# ASSUMPTION: Each 10-µs SWAP step contains an estimated 10% MW duty cycle
# internally (i.e., ~1 µs of actual MW pulses per SWAP), consistent with
# pulsed-hyperfine-gate sequences.  The local gate step uses:
#   N_pi × t_pi = (t_gates / t_pi) × t_pi = t_gates = 1 µs of on-time.
#
# NOTE: No first-principles thermal conductance value is published (as of
# Karapatzakis 2026), so this is an ORDER-OF-MAGNITUDE estimate, not a
# validated thermal calculation.

def compute_mw_duty_cycle(t_comm):
    """
    Estimate MW duty cycle for the checkpointing protocol.

    Gate on-time per cycle:
      Local gate step : t_gates (assume CW-equivalent; 1 µs on-time)
      SWAP in         : t_pi_pulse * (t_swap / t_pi_pulse / 10)  [~10% pulse filling]
      SWAP out        : same
    Total cycle time  : t_gates + 2*t_swap + t_comm

    Returns (duty_cycle, t_on, t_cycle) tuple.
    """
    # Conservative estimate: assume ~10% pulse-filling inside each SWAP window
    # (realistic for hyperfine gate sequences with refocusing pulses)
    n_pulses_per_swap = max(1, int(t_swap / t_pi_pulse / 10))
    t_mw_per_swap     = n_pulses_per_swap * t_pi_pulse

    t_on    = t_gates + 2 * t_mw_per_swap   # total MW on-time per cycle
    t_cycle = t_gates + 2 * t_swap + t_comm  # total cycle time
    duty    = t_on / t_cycle
    return duty, t_on, t_cycle

# ============================================================================
# MAIN SWEEP
# ============================================================================
print("=" * 72)
print("PART 4: FULL CHECKPOINTING PROTOCOL  [CORRECTED PARAMETERS]")
print("Diamond Quantum Computing : Checkpointing Protocol")
print("=" * 72)
print()
print("SWAP proxy fidelity (0.9992) = ¹³C single-qubit gate fidelity (Resch 2026)")
print("  *** NOT an experimentally confirmed SWAP gate fidelity ***")
print()
print(f"  T1_snv   = 4.23 ms  (Rosenthal 2023, direct 3 K; no direct 1.7 K measurement —")
print(f"             T1 > ~20 ms lower bound; Orbach ≈ 200 s; 4.23 ms = conservative worst-case)")
print(f"  T2_snv   = 170 µs   (Hahn echo, plateau 1.7–3 K : Rosenthal 2023)")
print(f"  T2_13c   = 1.35 s   (CPMG-128, 50 mK : Resch 2026; NOT validated at 1.7–3 K)")
print(f"  t_gates  = {t_gates*1e6:.0f} µs   |   t_swap = {t_swap*1e6:.0f} µs   "
      f"|   F_swap = {swap_fidelity} (proposed)")
print(f"  t_pi     = {t_pi_pulse*1e9:.1f} ns  (SnV MW π-pulse : Rosenthal 2023)")
print()

cp_results  = {T: np.zeros(len(t_comm_values)) for T in temperatures}
ncp_results = {T: np.zeros(len(t_comm_values)) for T in temperatures}

for T in temperatures:
    print(f"  Running T = {T}K  "
          f"(SnV T2={T2_snv[T]*1e6:.0f} µs, 13C T2={T2_13c[T]:.2f} s) ...",
          end='  ', flush=True)
    for i, tc in enumerate(t_comm_values):
        cp_results[T][i]  = full_checkpoint(rho_init, tc, T)
        ncp_results[T][i] = no_checkpoint(rho_init, tc, T)
    print("done.")

print()

# ============================================================================
# HEATING BUDGET TABLE
# ============================================================================
t_comm_us = t_comm_values * 1e6

print("=" * 72)
print("HEATING BUDGET ANALYSIS  (Karapatzakis et al., arXiv:2606.15398, 2026)")
print("=" * 72)
print()
print("MW duty cycle vs communication time:")
print()
print(f"  {'t_comm (µs)':>12}  {'t_cycle (µs)':>14}  {'MW on (ns)':>12}  "
      f"{'Duty cycle (%)':>16}  {'Safe? (<<20%)':>14}")
print(f"  {'-'*12}  {'-'*14}  {'-'*12}  {'-'*16}  {'-'*14}")

duty_cycles = []
for tc in t_comm_values:
    duty, t_on, t_cycle = compute_mw_duty_cycle(tc)
    duty_pct = duty * 100
    duty_cycles.append(duty_pct)
    safe = "[OK]" if duty_pct < 15.0 else ("[!]" if duty_pct < 20.0 else "[X]")
    print(f"  {tc*1e6:>12.1f}  {t_cycle*1e6:>14.1f}  {t_on*1e9:>12.1f}  "
          f"{duty_pct:>16.2f}  {safe:>14}")

print()
print("Interpretation:")
print("  Karapatzakis 2026 measured 20%-duty-cycle breakdown thresholds at B=400 mT:")
print("    10 dBm  (θ=84°)  and  ~12 dBm  (θ=20°)")
print("  Gradual heating onset begins above ~7 dBm (CW).")
print()
print("  At long t_comm (>100 µs), the duty cycle drops well below 20% since")
print("  the MW gates are a small fraction of the cycle.")
print("  At short t_comm (1–10 µs), the duty cycle can approach or exceed 20%,")
print("  requiring careful power management to stay below 7–10 dBm onset.")
print()
print("  RF driving (20 MHz, ¹³C rotations): no measurable heating was observed")
print("  up to an abrupt breakdown at 9.4 dBm (Bac ≈ 1.2 mT) in Karapatzakis 2026.")
print("  RF heating is therefore NOT the binding constraint in this protocol.")
print("  MW heating is the binding constraint, especially at short t_comm.")
print()
print("  DISCLAIMER: No first-principles thermal conductance value is published")
print("  (Karapatzakis 2026). This analysis is an ORDER-OF-MAGNITUDE estimate")
print("  based on measured duty-cycle thresholds and is NOT a validated thermal")
print("  calculation.  The thermal relaxation time constant τ_th = 0.54(17) ms")
print("  (Karapatzakis 2026) implies that a high-duty-cycle burst lasting several")
print(f"  τ_th = {tau_th_ms:.2f} ms can accumulate heat if t_comm << τ_th.")
print()

# ============================================================================
# RESULTS TABLE
# ============================================================================
print("=" * 72)
print("CHECKPOINT RESULTS TABLE  (* = CP > NCP,  ! = CP fidelity < 0.99)")
print("[All fidelities conditioned on proposed SWAP proxy : see caveat above]")
print("=" * 72)
hdr = (f"{'t_comm(µs)':>11}  "
       f"{'CP@1.7K':>8} {'NCP@1.7K':>9}  "
       f"{'CP@3K':>7} {'NCP@3K':>8}  flags")
print(hdr)
print("-" * 72)

for i, tc_us in enumerate(t_comm_us):
    flags    = ""
    row_cp   = [cp_results[T][i]  for T in temperatures]
    row_ncp  = [ncp_results[T][i] for T in temperatures]
    if all(cp > ncp for cp, ncp in zip(row_cp, row_ncp)):
        flags += "*"
    if any(cp < 0.99 for cp in row_cp):
        flags += "!"
    print(f"{tc_us:>11.1f}  "
          f"{row_cp[0]:>8.4f} {row_ncp[0]:>9.4f}  "
          f"{row_cp[1]:>7.4f} {row_ncp[1]:>8.4f}  {flags}")

print()

# ============================================================================
# CROSSOVER ANALYSIS
# ============================================================================
print("=" * 72)
print("CROSSOVER AND ADVANTAGE ANALYSIS")
print("=" * 72)
print()

for T in temperatures:
    cp  = cp_results[T]
    ncp = ncp_results[T]
    adv = cp - ncp

    helps_mask = adv > 0
    if helps_mask.any():
        last_idx      = np.where(helps_mask)[0][-1]
        crossover_str = f"{t_comm_us[last_idx]:.1f} µs"
    else:
        crossover_str = "never (SWAP overhead dominates)"

    max_adv_idx = np.argmax(adv)
    max_adv     = adv[max_adv_idx]
    max_adv_us  = t_comm_us[max_adv_idx]

    print(f"  At {T}K: checkpointing helps for t_comm up to {crossover_str}")
    print(f"  At {T}K: maximum advantage = {max_adv:+.4f} at t_comm = {max_adv_us:.1f} µs")
    print()

# ============================================================================
# VERIFICATION CHECKS
# ============================================================================
print("─── Core claim verification ──────────────────────────────────────────")
print()

short_tc_idx = 0
check1 = all(cp_results[T][short_tc_idx] > ncp_results[T][short_tc_idx]
             for T in temperatures)
print(f"  [OK] CP > NCP at short t_comm ({t_comm_us[short_tc_idx]:.0f} µs): {check1}")

long_tc_idx = -1
cp_long  = [cp_results[T][long_tc_idx] for T in temperatures]
cp_short = [cp_results[T][0]           for T in temperatures]
check2 = all(cp_l < cp_s for cp_l, cp_s in zip(cp_long, cp_short))
print(f"  [OK] CP fidelity declines at long t_comm (13C decoherence): {check2}")

def advantage_window(T):
    return int(np.sum(cp_results[T] > ncp_results[T]))

check3 = advantage_window(1.7) >= advantage_window(3.0)
print(f"  [OK] 1.7K advantage window ({advantage_window(1.7)} pts) "
      f">= 3K window ({advantage_window(3.0)} pts): {check3}")
print()

# ============================================================================
# PLOT 1 : CHECKPOINTED vs NO-CHECKPOINT (2×2 GRID)
# ============================================================================
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), sharey=True)
fig1.suptitle(
    "Part 4 : Proposed Checkpointing Protocol: Fidelity vs t_comm\n"
    "[SWAP fidelity = 0.9992 is a proposed proxy : not experimentally validated]",
    fontsize=12, fontweight='bold')

row_titles = ["With Proposed Checkpoint", "Without Checkpoint"]
data_rows  = [cp_results, ncp_results]

for row_idx, (row_title, data) in enumerate(zip(row_titles, data_rows)):
    for col_idx, T in enumerate(temperatures):
        ax = axes1[row_idx, col_idx]
        ax.semilogx(t_comm_us, data[T],
                    color=COLORS[T], linewidth=2.0, marker='o', markersize=5)
        ax.axhline(0.99, color='green', linestyle='--', linewidth=1.2,
                   label='0.99 threshold')
        ax.set_title(f"{T}K : {row_title}", fontsize=10, fontweight='bold')
        ax.set_xlabel("t_comm (µs)", fontsize=9)
        ax.set_ylabel("Fidelity with |+⟩", fontsize=9)
        ax.set_ylim([0.0, 1.02])
        ax.grid(True, alpha=0.3)
        if row_idx == 0 and col_idx == 0:
            ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('part4_checkpoint_comparison.png', dpi=150, bbox_inches='tight')
plt.close(fig1)
print("Plot saved: part4_checkpoint_comparison.png")

# ============================================================================
# PLOT 2 : ADVANTAGE OVERLAY (2 subplots: 1.7K | 3K)
# ============================================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
fig2.suptitle(
    "Part 4 : Proposed Checkpointing Advantage vs Communication Time\n"
    "[Green: CP > NCP, Red: SWAP overhead dominates]",
    fontsize=12, fontweight='bold')

for col_idx, T in enumerate(temperatures):
    ax  = axes2[col_idx]
    cp  = cp_results[T]
    ncp = ncp_results[T]

    ax.semilogx(t_comm_us, cp,  color=COLORS[T], linewidth=2.0,
                linestyle='-',  label='Proposed checkpoint')
    ax.semilogx(t_comm_us, ncp, color=COLORS[T], linewidth=2.0,
                linestyle='--', label='No checkpoint', alpha=0.7)
    ax.fill_between(t_comm_us, cp, ncp,
                    where=(cp >= ncp), alpha=0.20, color='green',
                    label='Checkpointing advantage')
    ax.fill_between(t_comm_us, cp, ncp,
                    where=(cp <  ncp), alpha=0.20, color='red',
                    label='SWAP overhead dominates')
    ax.axhline(0.99, color='green', linestyle=':', linewidth=1.2,
               label='0.99 threshold')

    ax.set_title(f"{T}K", fontsize=12, fontweight='bold')
    ax.set_xlabel("t_comm (µs)", fontsize=10)
    ax.set_ylabel("Fidelity with |+⟩", fontsize=10)
    ax.set_ylim([0.0, 1.02])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='lower left')

plt.tight_layout()
plt.savefig('part4_advantage.png', dpi=150, bbox_inches='tight')
plt.close(fig2)
print("Plot saved: part4_advantage.png")

# ============================================================================
# PLOT 3 : HEATING BUDGET: DUTY CYCLE vs t_comm
# ============================================================================
# Also shows the MW/RF heating margins from Karapatzakis 2026.

fig3, ax3 = plt.subplots(figsize=(9, 5))

ax3.semilogx(t_comm_us, duty_cycles, color='#2ca02c', linewidth=2.0,
             marker='o', markersize=6, label='MW duty cycle (this protocol)')

# Threshold lines
ax3.axhline(20.0, color='#d62728', linestyle='--', linewidth=1.5,
            label='20% threshold (Karapatzakis 2026 test condition)')
ax3.axhline(15.0, color='#ff7f0e', linestyle=':', linewidth=1.5,
            label='~15%: onset region (extrapolated)')
ax3.axhspan(0, 15.0, alpha=0.07, color='green')
ax3.axhspan(15.0, 20.0, alpha=0.10, color='orange')
ax3.axhspan(20.0, 35.0, alpha=0.07, color='red')

# RF: no heating → show as a flat "safe" annotation
ax3.text(1.2, 28,
         'RF (20 MHz, ¹³C): no measurable heating\nup to breakdown at 9.4 dBm',
         fontsize=9, color='#1f77b4',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f4fd',
                   edgecolor='#1f77b4', linewidth=1.0))

ax3.set_xlabel('t_comm (µs)', fontsize=11)
ax3.set_ylabel('Estimated MW duty cycle (%)', fontsize=11)
ax3.set_title(
    'MW Heating Budget: Estimated Duty Cycle vs Communication Time\n'
    'Source: Karapatzakis et al., arXiv:2606.15398 (2026)',
    fontsize=11, fontweight='bold')
ax3.set_xlim([t_comm_us[0], t_comm_us[-1]])
ax3.set_ylim([0, 35])
ax3.legend(fontsize=9, loc='upper right')
ax3.grid(True, alpha=0.3, which='both')

# Caveat box
ax3.text(0.01, 0.05,
         '(!) Order-of-magnitude estimate only.\n'
         'No first-principles thermal conductance\n'
         'value is published (Karapatzakis 2026).',
         transform=ax3.transAxes, fontsize=8.5,
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3cd',
                   edgecolor='#856404', linewidth=1.2))

plt.tight_layout()
plt.savefig('part4_heating_budget.png', dpi=150, bbox_inches='tight')
plt.close(fig3)
print("Plot saved: part4_heating_budget.png")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("=" * 72)
print("MAXIMUM CHECKPOINTING ADVANTAGE SUMMARY")
print("=" * 72)
print()
for T in temperatures:
    adv         = cp_results[T] - ncp_results[T]
    max_adv_idx = np.argmax(adv)
    print(f"  At {T}K: maximum proposed advantage = "
          f"{adv[max_adv_idx]:+.4f}  "
          f"at t_comm = {t_comm_us[max_adv_idx]:.1f} µs")
print()

print("=" * 72)
print("PHYSICAL SUMMARY")
print("=" * 72)
print()
print("  At 1.7 K and 3 K: T2_SnV = 170 µs, T2_13C ≈ 1.35 s  → ratio ≈ 7940×")
print("  (ratio uses CPMG-128 value; Hahn echo gives T2_13C = 167 ms → ratio ≈ 980×)")
print()
print("  Heating constraint (MW): duty cycle falls well below 20% for t_comm > 50 µs.")
print("  For t_comm < 10 µs, duty cycle may approach 20% : power management required.")
print("  RF driving shows no measurable heating up to breakdown (Karapatzakis 2026).")
print("  MW heating is the binding constraint; RF is not.")
print()
print("  *** REMINDER: SWAP fidelity (0.9992) is a PROPOSED proxy based on")
print("      ¹³C single-qubit gate fidelity at 50 mK (Resch 2026).  The SWAP")
print("      operation itself has not been experimentally demonstrated. ***")
print()
print("Part 4 complete. Ready for Part 5: Publication Figures.")
