"""
Part 2: Decoherence Simulation
===============================
Diamond quantum computing architecture : quantum checkpointing protocol.

Qubit 1: SnV electron spin  (compute qubit)
Qubit 2: 13C nuclear spin   (memory qubit)

This file covers:
  - Lindblad collapse operators for T1 (relaxation) and T2 (dephasing)
  - Lindblad master equation via qutip.mesolve()
  - Free evolution (H = 0) at 1.7 K and 3 K  [7 K and 10 K REMOVED —
    not supported by the published SnV data (see parameter notes below)]
  - Fidelity decay curves for both qubits

PHYSICAL PARAMETER SOURCES
----------------------------
SnV electron spin : Rosenthal et al., Phys. Rev. X 13, 031022 (2023)
  T1   = 4.23 ± 1.37 ms  measured at 3 K (Rosenthal 2023).
         No direct 1.7 K measurement exists; Rosenthal 2023 gives only a T1 > ~20 ms
         lower bound at 1.7 K, and the Orbach extrapolation predicts T1 ≈ 200 s there.
         T1_snv[1.7] = 4.23e-3 is retained as a conservative worst-case bound.
  NOTE: T1 = 5.22 ± 1.54 µs at 5 K : used only if 5 K simulations are requested
  T2_echo = 170.0 ± 2.8 µs  (Hahn echo, plateau from 1.7 K to ~3 K)
         Extendable to 650 µs via XY16 dynamical decoupling (not used here)
  T2*     = 396.6 ± 2.29 ns  (Ramsey, Fig. 3e at 1.7 K; inhomogeneous/nuclear-bath-limited)
  Orbach: Γph(T) ∝ Δg³ / (exp(α·ℏΔg / kT) − 1),  Δg/2π = 903.0 GHz, α = 1.207

13C nuclear spin : Resch et al., Phys. Rev. X 16, 011060 (2026) [arXiv:2509.03354]
  T2*        = 1.5(1) ms    (Ramsey)
  T2_Hahn    = 167(9) ms   (Hahn echo)
  T2_CPMG128 = 1.35(3) s   (128-pulse CPMG : used as memory coherence time here)
  Init. fidelity: 99.74(3)%,  RF gate fidelity: 99.92(1)%

  *** IMPORTANT CAVEAT ***
  These 13C values were measured at 50 mK in a dilution refrigerator.
  They have NOT been demonstrated at 1.7–3 K.  No SWAP gate was demonstrated
  in this source : initialization used optical+MW pumping, not a SWAP operation.
  The SWAP fidelity used in Parts 3–4 (0.9992) is the single-qubit RF gate
  fidelity from Resch et al., used here as a proposed/unvalidated SWAP proxy.

No SWAP. No checkpointing. No driving Hamiltonian.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, qeye, sigmaz, sigmap, sigmam,
    ket2dm, ptrace, mesolve, expect
)

# ============================================================================
# SYSTEM REBUILD  (standalone : no import from Part 1)
# ============================================================================

I2 = qeye(2)

snv_0   = basis(2, 0)
snv_1   = basis(2, 1)
psi_snv = (snv_0 + snv_1).unit()   # |+>_SnV

c13_0   = basis(2, 0)
psi_13c = c13_0                    # |0>_13C

psi_full = tensor(psi_snv, psi_13c)
rho0     = ket2dm(psi_full)        # 4×4 density matrix at t=0

# ============================================================================
# PHYSICAL PARAMETERS  : Corrected, literature-traceable values
# ============================================================================
# All times in seconds.
#
# SnV T1 at operating temperatures (Rosenthal et al., PRX 13, 031022 (2023))
# The Orbach fit gives a flat plateau in T1 from ~1.7 K to just below 3 K.
# Both temperatures share T1 = 4.23 ms to within measurement uncertainty.
T1_snv = {
    1.7: 4.23e-3,   # 4.23 ms  : Rosenthal 2023, 3 K value used as conservative worst-case
                    #             (no direct 1.7 K measurement; T1 > ~20 ms lower bound;
                    #              Orbach extrapolation ≈ 200 s at 1.7 K)
    3.0: 4.23e-3,   # 4.23 ms  : Rosenthal 2023 direct measurement at 3 K
}

# SnV T2 (Hahn echo, plateau) : Rosenthal et al., PRX 13, 031022 (2023)
# T2_echo = 170 µs is measured as flat (within error) across 1.7 K to ~3 K.
# T2* = 396.6 ± 2.29 ns (Rosenthal 2023, Fig. 3e Ramsey fit at 1.7 K) is NOT used for
# the Lindblad simulation (it is the inhomogeneous free-induction-decay time and would
# require a quasi-static bath model); T2_echo is the physically meaningful input.
T2_snv = {
    1.7: 170.0e-6,   # 170 µs  : Rosenthal 2023
    3.0: 170.0e-6,   # 170 µs  : Rosenthal 2023 (same plateau)
}

# 13C T2 : Resch et al., PRX 16, 011060 (2026)
# We use the CPMG-128 value (1.35 s) as the memory coherence input.
# This is the most relevant figure of merit for a storage qubit.
# *** Measured at 50 mK : not validated at 1.7–3 K (see module docstring). ***
T2_13c = {
    1.7: 1.35,   # 1.35 s  : Resch 2026, 128-pulse CPMG at 50 mK
    3.0: 1.35,   # same value; no temperature dependence data exists at 1.7–3 K
}

# Simulation temperatures (Kelvin) : ONLY 1.7 K and 3 K supported by data
temperatures = [1.7, 3.0]

# Colour scheme for plots
COLORS = {1.7: '#1f77b4', 3.0: '#d62728'}   # blue = 1.7 K, red = 3 K

# ============================================================================
# OPERATORS ON THE FULL HILBERT SPACE
# ============================================================================
Sz_snv = tensor(sigmaz(), I2)
Sz_13c = tensor(I2, sigmaz())
Sm_snv = tensor(sigmam(), I2)

# ============================================================================
# TIME ARRAY
# ============================================================================
# Simulate 0 → 500 µs  (3× T2_snv = 3×170 µs = 510 µs, well within window).
# At 3 ms the SnV decoherence is complete but the window would be dominated
# by the flat 13C curve; 500 µs gives more informative SnV detail.
t_start = 0.0
t_end   = 500e-6    # 500 µs
n_steps = 300
times   = np.linspace(t_start, t_end, n_steps)

# ============================================================================
# HAMILTONIAN
# ============================================================================
H = tensor(qeye(2), qeye(2)) * 0.0

# ============================================================================
# LINDBLAD COLLAPSE OPERATORS
# ============================================================================
# Correct prefactor convention (verified, unchanged from prior version):
#
#   For C = sqrt(gamma) * sigmaz the Lindblad dissipator gives:
#     d/dt rho_01 = -2*gamma * rho_01
#   so the physical dephasing time is T2 = 1/(2*gamma).
#   To target T2_echo: need gamma_phi = 1/(2*T2_echo).
#
#   The amplitude-damping operator Sm contributes additional dephasing
#   at rate 1/(2*T1).  The TOTAL coherence decay rate is:
#     1/T2_eff = 2*gamma_phi + 1/(2*T1)
#   Experimental T2 = 1/( gamma_phi + 1/(2*T1) ) ... so:
#     gamma_phi = 1/T2 - 1/(2*T1)   (pure-dephasing contribution only)
#
#   If T2 > 2*T1 (violates quantum bound T2 ≤ 2*T1), set gamma_phi = 0
#   (system is T1-limited; T2_eff = 2*T1).
#   At 1.7 K: T2 = 170 µs, 2*T1 = 8.46 ms → T2 < 2*T1 ✓ (not T1-limited)
#   At 3.0 K: T2 = 170 µs, 2*T1 = 8.46 ms → T2 < 2*T1 ✓ (not T1-limited)
#
#   13C has no measurable T1 over this window, so its operator is
#   sqrt(1/(2*T2_13c)) * Sz_13c only.

# ============================================================================
# LINDBLAD MASTER EQUATION SIMULATION
# ============================================================================
fid_snv_vs_time = {}
fid_13c_vs_time = {}

print("=" * 70)
print("PART 2: DECOHERENCE SIMULATION  [CORRECTED PARAMETERS]")
print("Diamond Quantum Computing : Checkpointing Protocol")
print("=" * 70)
print()
print("Physical parameters (literature-traceable):")
print("  SnV T1    : 4.23 ms at 1.7 K and 3 K  (Rosenthal et al., PRX 2023)")
print("  SnV T2    : 170 µs  (Hahn echo, plateau 1.7–3 K) (Rosenthal 2023)")
print("  13C T2    : 1.35 s  (128-pulse CPMG at 50 mK) (Resch et al., PRX 2026)")
print("  *** CAVEAT: 13C data at 50 mK : not validated at 1.7–3 K ***")
print()
print(f"Simulation window : 0 – {t_end*1e6:.0f} µs  ({n_steps} steps)")
print(f"Initial SnV state : |+>  (Bloch vector along +x)")
print(f"Initial 13C state : |0>  (Bloch vector along +z)")
print()

for T in temperatures:
    print(f"Running mesolve() at T = {T}K  "
          f"(SnV T1={T1_snv[T]*1e3:.2f} ms, SnV T2={T2_snv[T]*1e6:.0f} µs, "
          f"13C T2={T2_13c[T]:.2f} s) ...",
          end='  ', flush=True)

    gamma_phi_snv = max(0.0, 1.0/T2_snv[T] - 1.0/(2*T1_snv[T]))
    c_ops = [
        np.sqrt(gamma_phi_snv / 2.0)     * Sz_snv,   # SnV pure dephasing
        np.sqrt(1.0 / T1_snv[T])         * Sm_snv,   # SnV relaxation
        np.sqrt(1.0 / (2*T2_13c[T]))     * Sz_13c,   # 13C dephasing
    ]

    result = mesolve(H, rho0, times, c_ops, e_ops=[])

    fid_snv_t = np.zeros(n_steps)
    fid_13c_t = np.zeros(n_steps)

    for i, rho_t in enumerate(result.states):
        rho_snv_t = ptrace(rho_t, 0)
        rho_13c_t = ptrace(rho_t, 1)
        fid_snv_t[i] = float((psi_snv.dag() * rho_snv_t * psi_snv).real)
        fid_13c_t[i] = float((psi_13c.dag() * rho_13c_t * psi_13c).real)

    fid_snv_vs_time[T] = fid_snv_t
    fid_13c_vs_time[T] = fid_13c_t
    print("done.")

print()

# ============================================================================
# SUMMARY TABLE
# ============================================================================
# Time points in microseconds
table_times_us = [0, 50, 100, 170, 250, 350, 500]

def time_index(t_us):
    t_s = t_us * 1e-6
    return int(np.argmin(np.abs(times - t_s)))

print("=" * 70)
print("SUMMARY TABLE: SnV FIDELITY vs TIME")
print("=" * 70)
header = f"{'Time (µs)':>10}  {'SnV F@1.7K':>12}  {'SnV F@3K':>10}"
print(header)
print("-" * 42)
for t_us in table_times_us:
    idx = time_index(t_us)
    row = f"{t_us:>10}"
    for T in temperatures:
        row += f"  {fid_snv_vs_time[T][idx]:>12.6f}"
    print(row)
print()

# ============================================================================
# INTERPRETATION: T2 CROSSING TIMES
# ============================================================================
print("=" * 70)
print("INTERPRETATION")
print("=" * 70)
print()

for T in temperatures:
    fid_arr  = fid_snv_vs_time[T]
    times_us = times * 1e6

    idx_99 = np.argmax(fid_arr < 0.99)
    t_99_str = (f"{times_us[idx_99]:.1f} µs"
                if fid_arr[idx_99] < 0.99
                else f"> {t_end*1e6:.0f} µs (within simulation window)")

    idx_50 = np.argmax(fid_arr < 0.50)
    t_50_str = (f"{times_us[idx_50]:.1f} µs"
                if fid_arr[idx_50] < 0.50
                else f"> {t_end*1e6:.0f} µs (within simulation window)")

    print(f"  At {T}K: SnV reaches 99% fidelity threshold at t = {t_99_str}")
    print(f"  At {T}K: SnV fully decoheres (F<0.50)   at t = {t_50_str}")
    print()

# ============================================================================
# PLOT 1 : SnV FIDELITY vs TIME
# ============================================================================
times_us = times * 1e6   # µs for x-axis

fig1, ax1 = plt.subplots(figsize=(10, 6))

T2_labels = {1.7: 'T₂ = 170 µs', 3.0: 'T₂ = 170 µs'}
# 1.7 K and 3 K are indistinguishable (same T1, T2) : plot a single plateau curve.
fid_arr = fid_snv_vs_time[1.7]
ax1.plot(times_us, fid_arr,
         color='#1f77b4', linewidth=2.0,
         label='1.7–3 K  (T₂ = 170 µs, Rosenthal 2023)\nCurves coincide : no resolvable difference in this range')

# Mark single T2 point
t2_us = T2_snv[1.7] * 1e6
if t2_us <= times_us[-1]:
    idx_t2 = int(np.argmin(np.abs(times - T2_snv[1.7])))
    ax1.plot(times_us[idx_t2], fid_arr[idx_t2],
             'o', color='#1f77b4', markersize=9, zorder=5,
             label=f'T₂ = {t2_us:.0f} µs (1.7–3 K plateau)')

# Reference lines
ax1.axhline(0.99, color='green',  linestyle='--', linewidth=1.5,
            label='Fault tolerance threshold (0.99)')
ax1.axhline(0.50, color='grey',   linestyle='--', linewidth=1.5,
            label='Fully decohered (0.50)')

ax1.set_xlabel('Time (µs)', fontsize=12)
ax1.set_ylabel('Fidelity with |+⟩', fontsize=12)
ax1.set_title(
    'SnV Electron Spin Decoherence : 1.7–3 K Plateau\n'
    '(T₂ = 170 µs Hahn echo : Rosenthal et al., PRX 2023)',
    fontsize=12, fontweight='bold')
ax1.set_xlim([0, t_end * 1e6])
ax1.set_ylim([0.0, 1.02])
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('part2_snv_decoherence.png', dpi=150, bbox_inches='tight')
plt.close(fig1)
print("Plot saved: part2_snv_decoherence.png")

# ============================================================================
# PLOT 2 : 13C FIDELITY vs TIME
# ============================================================================
fig2, ax2 = plt.subplots(figsize=(10, 6))

# 1.7 K and 3 K are indistinguishable : plot a single plateau curve.
ax2.plot(times_us, fid_13c_vs_time[1.7],
         color='#1f77b4', linewidth=2.0, linestyle='--',
         label='1.7–3 K  (T₂(CPMG-128) = 1.35 s : Resch 2026, 50 mK)\nCurves coincide : no resolvable difference in this range')

# Reference lines
ax2.axhline(0.99, color='green', linestyle='--', linewidth=1.5,
            label='Fault tolerance threshold (0.99)')
ax2.axhline(0.50, color='grey',  linestyle='--', linewidth=1.5,
            label='Fully decohered (0.50)')

ax2.set_xlabel('Time (µs)', fontsize=12)
ax2.set_ylabel('Fidelity with |0⟩', fontsize=12)
ax2.set_title(
    '¹³C Nuclear Spin Decoherence : 1.7–3 K Plateau\n'
    '(T₂(CPMG-128) = 1.35 s : Resch et al., PRX 2026 @ 50 mK)',
    fontsize=12, fontweight='bold')
ax2.set_xlim([0, t_end * 1e6])
ax2.set_ylim([0.9998, 1.00002])   # 13C barely moves over 500 µs
ax2.legend(fontsize=9, loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.text(
    0.98, 0.05,
    '(!) CAVEAT: 13C values measured at 50 mK.\n'
    '   Not validated at 1.7-3 K operating temperature.',
    transform=ax2.transAxes, fontsize=9, ha='right', va='bottom',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3cd',
              edgecolor='#856404', linewidth=1.2))
plt.tight_layout()
plt.savefig('part2_13c_decoherence.png', dpi=150, bbox_inches='tight')
plt.close(fig2)
print("Plot saved: part2_13c_decoherence.png")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("=" * 70)
print("PHYSICAL INTERPRETATION SUMMARY")
print("=" * 70)
print()
print("SnV electron spin (Rosenthal et al., PRX 2023):")
print("  - T2 = 170 µs (Hahn echo) across 1.7–3 K : nearly temperature-independent")
print("  - T2* = 396.6 ± 2.29 ns (Rosenthal 2023, Fig. 3e Ramsey fit at 1.7 K; not used in Lindblad model)")
print("  - T1 = 4.23 ms used at both temperatures (Rosenthal 2023, direct 3 K measurement)")
print("       No direct 1.7 K measurement; T1 > ~20 ms lower bound; Orbach ≈ 200 s.")
print("       4.23 ms at 1.7 K is a conservative worst-case bound.")
print("       This choice does not affect simulation results: T2_echo = 170 µs")
print("       dominates coherence decay at both temperatures (T1 >> T2_echo).")
print("  - Dephasing is dominated by pure dephasing (gamma_phi >> 1/(2*T1))")
print()
print("13C nuclear spin (Resch et al., PRX 2026, at 50 mK):")
print("  - T2(CPMG-128) = 1.35 s : orders of magnitude longer than SnV T2")
print("  - Over the 500 µs window, 13C fidelity is essentially 1.000")
print("  - *** NOT validated at 1.7–3 K. SWAP protocol not experimentally confirmed. ***")
print()
print("Conclusion: the enormous T2_13C / T2_SnV ≈ 8000 ratio (at CPMG-128 values)")
print("  motivates the checkpointing protocol : even with significant caveats.")
print()
print("Part 2 complete. Ready for Part 3: SWAP gate.")
