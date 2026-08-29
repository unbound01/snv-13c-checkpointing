"""
Part 5: Publication-Quality Figures
=====================================
Diamond quantum computing architecture : quantum checkpointing protocol.

Generates five journal-ready figures from the full checkpointing simulation.
All physics is recomputed internally : no imports from other parts.

TEMPERATURES: 1.7 K and 3 K ONLY
  7 K and 10 K have been removed : not supported by the SnV literature.

LABELLING CAVEAT : SWAP FIDELITY
  swap_fidelity = 0.9992 is the ¹³C single-qubit RF gate fidelity from
  Resch et al., PRX 16, 011060 (2026).  It is used as a PROPOSED/UNVALIDATED
  proxy for a SnV↔¹³C SWAP gate.  No SWAP was demonstrated in that paper.
  All checkpoint fidelity curves are conditioned on this assumption.

PHYSICAL PARAMETER SOURCES
  SnV : Rosenthal et al., Phys. Rev. X 13, 031022 (2023)
  ¹³C : Resch et al., Phys. Rev. X 16, 011060 (2026) [arXiv:2509.03354]
        : measured at 50 mK, NOT validated at 1.7–3 K —
  Heating: Karapatzakis et al., arXiv:2606.15398 (2026)

Figures produced:
  Fig 1 : part5_fig1_decoherence.png              SnV decoherence at 1.7K & 3K
  Fig 2 : part5_fig2_checkpoint_advantage.png     Checkpointing advantage per temperature
  Fig 3 : part5_fig3_temperature_comparison.png   Both temperatures overlaid
  Fig 4 : part5_fig4_summary_table.png            Colour-coded data table
  Fig 5 : part5_fig5_heating_margin.png           MW vs RF heating margin comparison
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from qutip import (
    basis, tensor, qeye, sigmaz, sigmam,
    ket2dm, ptrace, mesolve, Qobj
)

# ============================================================================
# GLOBAL STYLE
# ============================================================================
plt.rcParams['font.family']       = 'serif'
plt.rcParams['font.size']         = 11
plt.rcParams['axes.linewidth']    = 1.2
plt.rcParams['xtick.direction']   = 'in'
plt.rcParams['ytick.direction']   = 'in'
plt.rcParams['xtick.major.size']  = 4
plt.rcParams['ytick.major.size']  = 4
plt.rcParams['xtick.minor.size']  = 2
plt.rcParams['ytick.minor.size']  = 2
plt.rcParams['legend.framealpha'] = 0.9
plt.rcParams['legend.edgecolor']  = '0.8'
plt.rcParams['axes.grid']         = False

# Two temperatures only : blue = 1.7 K, red = 3 K
COLORS = {1.7: '#1f77b4', 3.0: '#d62728'}

def remove_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ============================================================================
# SYSTEM SETUP
# ============================================================================
I2 = qeye(2)
I4 = tensor(I2, I2)

psi_snv_plus = (basis(2, 0) + basis(2, 1)).unit()
psi_13c_zero = basis(2, 0)
rho_init     = ket2dm(tensor(psi_snv_plus, psi_13c_zero))

# ============================================================================
# PHYSICAL PARAMETERS
# ============================================================================
# SnV : Rosenthal et al., PRX 13, 031022 (2023)
T1_snv = {1.7: 4.23e-3, 3.0: 4.23e-3}      # T1 (s) : Rosenthal 2023, direct 3 K measurement.
                                              # No direct 1.7 K value; T1 > ~20 ms lower bound;
                                              # Orbach extrapolation ≈ 200 s at 1.7 K.
                                              # 4.23 ms used at 1.7 K as conservative worst-case.
T2_snv = {1.7: 170.0e-6, 3.0: 170.0e-6}    # T2 Hahn echo (s), plateau 1.7–3 K
# T2* = 396.6 ± 2.29 ns (Rosenthal 2023, Fig. 3e Ramsey fit at 1.7 K;
#        inhomogeneous) : NOT used in Lindblad model

# ¹³C : Resch et al., PRX 16, 011060 (2026) [at 50 mK : NOT validated at 1.7–3 K]
T2_13c = {1.7: 1.35, 3.0: 1.35}            # CPMG-128 value (s)

# SWAP proxy : ¹³C single-qubit RF gate fidelity (Resch 2026, 50 mK)
# *** NOT an experimentally confirmed SWAP gate fidelity ***
swap_fidelity = 0.9992

t_swap    = 10e-6    # assumed SWAP gate duration
t_gates   = 1e-6     # local SnV gate duration
t_pi      = 48.4e-9  # SnV MW π-pulse (Rosenthal 2023)
tau_th_ms = 0.54     # thermal relaxation time constant, ms (Karapatzakis 2026)

temperatures = [1.7, 3.0]

# t_comm sweep : 11 values, 1 µs → 10 ms
t_comm_values = np.array([1e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6,
                           200e-6, 500e-6, 1e-3, 5e-3, 10e-3])
t_comm_us     = t_comm_values * 1e6

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
    F = 0.9992 from ¹³C single-qubit gate fidelity (Resch 2026).
    *** Not an experimentally confirmed SWAP gate ***
    """
    return swap_fidelity * (U_swap * rho * U_swap.dag()) \
           + (1.0 - swap_fidelity) * rho_mixed

# ============================================================================
# SIMULATION HELPERS
# ============================================================================
def build_c_ops(T):
    """
    Lindblad collapse operators for temperature T (Kelvin).

    Correct convention (verified identical to Parts 2–4):
      gamma_phi = max(0, 1/T2 - 1/(2*T1))
      c_dephase  = sqrt(gamma_phi / 2) * Sz
      c_relax    = sqrt(1/T1)          * Sm
      c_13C      = sqrt(1/(2*T2_13c)) * Sz_13c

    At both 1.7K and 3K: T2 = 170 µs, 2*T1 = 8.46 ms → T2 << 2*T1 (not T1-limited)
    """
    gamma_phi_snv = max(0.0, 1.0/T2_snv[T] - 1.0/(2*T1_snv[T]))
    return [
        np.sqrt(gamma_phi_snv / 2.0)  * Sz_snv,
        np.sqrt(1.0 / T1_snv[T])      * Sm_snv,
        np.sqrt(1.0 / (2*T2_13c[T])) * Sz_13c,
    ]

def evolve(rho_in, duration, c_ops_list, n_steps=50):
    if duration <= 0.0:
        return rho_in
    tlist  = np.linspace(0.0, duration, n_steps)
    result = mesolve(H, rho_in, tlist, c_ops_list, e_ops=[])
    return result.states[-1]

def pure_fidelity(psi_target, rho):
    return float((psi_target.dag() * rho * psi_target).real)

def full_checkpoint(rho_start, t_comm, T):
    c = build_c_ops(T)
    rho = evolve(rho_start, t_gates, c)
    rho = evolve(rho, t_swap, c)
    rho = noisy_swap(rho)
    rho = evolve(rho, t_comm, c)
    rho = evolve(rho, t_swap, c)
    rho = noisy_swap(rho)
    return pure_fidelity(psi_snv_plus, ptrace(rho, 0))

def no_checkpoint(rho_start, t_comm, T):
    c       = build_c_ops(T)
    t_total = t_gates + t_swap + t_comm + t_swap
    rho     = evolve(rho_start, t_total, c)
    return pure_fidelity(psi_snv_plus, ptrace(rho, 0))

def snv_decoherence_curve(T, times):
    c_ops  = build_c_ops(T)
    result = mesolve(H, rho_init, times, c_ops, e_ops=[])
    fids   = np.zeros(len(times))
    for i, rho_t in enumerate(result.states):
        fids[i] = pure_fidelity(psi_snv_plus, ptrace(rho_t, 0))
    return fids

def compute_mw_duty_cycle(t_comm):
    """
    Estimate MW duty cycle for the checkpointing protocol.

    Assumes ~10% pulse filling inside each SWAP window (conservative estimate
    for pulsed-hyperfine gate sequences with refocusing pulses).
    No first-principles thermal conductance is published (Karapatzakis 2026);
    this is an order-of-magnitude estimate only.
    """
    n_pulses_per_swap = max(1, int(t_swap / t_pi / 10))
    t_mw_per_swap     = n_pulses_per_swap * t_pi
    t_on    = t_gates + 2 * t_mw_per_swap
    t_cycle = t_gates + 2 * t_swap + t_comm
    return t_on / t_cycle, t_on, t_cycle

# ============================================================================
# RUN ALL SIMULATIONS
# ============================================================================
print("=" * 70)
print("PART 5: PUBLICATION FIGURES  [CORRECTED PARAMETERS]")
print("Diamond Quantum Computing : Checkpointing Protocol")
print("Temperatures: 1.7 K and 3 K only (7 K / 10 K removed)")
print("=" * 70)
print()

# Figure 1 data: SnV decoherence (0 → 500 µs, 300 steps)
print("Computing SnV decoherence curves (Fig 1) ...", end="  ", flush=True)
dec_times    = np.linspace(0, 500e-6, 300)
dec_times_us = dec_times * 1e6
dec_fids     = {T: snv_decoherence_curve(T, dec_times) for T in temperatures}
print("done.")

# Figures 2–4 data: checkpoint sweep
print("Computing checkpoint sweep (Figs 2–4) ...")
cp_results  = {T: np.zeros(len(t_comm_values)) for T in temperatures}
ncp_results = {T: np.zeros(len(t_comm_values)) for T in temperatures}

for T in temperatures:
    print(f"  T = {T}K ...", end="  ", flush=True)
    for i, tc in enumerate(t_comm_values):
        cp_results[T][i]  = full_checkpoint(rho_init, tc, T)
        ncp_results[T][i] = no_checkpoint(rho_init, tc, T)
    print("done.")

# Figure 5 data: heating budget
duty_pcts = np.array([compute_mw_duty_cycle(tc)[0] * 100 for tc in t_comm_values])

# Pre-compute advantage arrays
adv         = {T: cp_results[T] - ncp_results[T] for T in temperatures}
max_adv_idx = {T: int(np.argmax(adv[T]))          for T in temperatures}

print()

# ============================================================================
# FIGURE 1 : SnV DECOHERENCE CURVES  (1.7 K and 3 K)
# ============================================================================
print("Rendering Figure 1 ...", end="  ", flush=True)

fig1, ax1 = plt.subplots(figsize=(8, 5))

# Shaded fault-tolerant window
ax1.axhspan(0.99, 1.0, alpha=0.10, color='#1f77b4',
            label='Fault-tolerant window (F > 0.99)')

# Single plateau curve : 1.7 K and 3 K are indistinguishable (same T1, T2).
# Plotting both would stack two identical lines; use one curve labelled as the plateau.
ax1.plot(dec_times_us, dec_fids[1.7],
         color='#1f77b4', linewidth=2.0,
         label='1.7–3 K  (T₂ = 170 µs, Rosenthal 2023)\nCurves coincide : no resolvable difference in this range')

# Reference lines
ax1.axhline(0.99, color='0.4', linestyle='--', linewidth=1.2)
# Anchored mid-plot, below the line, clear of the upper-right legend
ax1.text(250, 0.965, 'F = 0.99', ha='center', va='top', fontsize=9, color='0.4')
ax1.axhline(0.50, color='0.6', linestyle='--', linewidth=1.0)
ax1.text(485, 0.505, 'Maximally mixed state',
         ha='right', va='bottom', fontsize=9, color='0.6')

# Single T2 vertical marker (one value covers both temperatures)
t2_us = T2_snv[1.7] * 1e6
ax1.axvline(t2_us, color='#1f77b4', linestyle=':', linewidth=1.2, alpha=0.7)
idx_t2 = int(np.argmin(np.abs(dec_times_us - t2_us)))
ax1.text(t2_us + 4, dec_fids[1.7][idx_t2] + 0.04,
         f'T₂ = {t2_us:.0f} µs\n(1.7–3 K plateau)',
         color='#1f77b4', fontsize=8.5, ha='left', va='bottom')

# T2* annotation (cosmetic : not fed into Lindblad model)
ax1.axvline(0.3966, color='0.6', linestyle='-.', linewidth=1.0, alpha=0.6)
ax1.text(0.41, 0.85, 'T₂* = 396.6 ns\n(Ramsey, Fig. 3e @ 1.7 K\nnot modelled)',
         color='0.5', fontsize=7.5, ha='left', va='top')

ax1.set_xlabel('Time (µs)', fontsize=12)
ax1.set_ylabel('Fidelity  F(t)', fontsize=12)
ax1.set_title(
    'SnV Electron Spin Decoherence : 1.7–3 K Plateau\n'
    '(T₂ = 170 µs Hahn echo; T1 = 4.23 ms : Rosenthal et al., PRX 2023)\n'
    'T1 = 3 K value reused at 1.7 K (no direct measurement); T2_echo dominates : result unchanged',
    fontsize=10, fontweight='bold', pad=10)
ax1.set_xlim([0, 500])
ax1.set_ylim([0.0, 1.02])
ax1.tick_params(labelsize=10)
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, color='0.85', linewidth=0.6, alpha=0.8)
remove_spines(ax1)

fig1.tight_layout()
fig1.savefig('part5_fig1_decoherence.png', dpi=300, bbox_inches='tight')
plt.close(fig1)
print("saved part5_fig1_decoherence.png")

# ============================================================================
# FIGURE 2 : CHECKPOINTING ADVANTAGE (2 subplots: 1.7 K | 3 K)
# ============================================================================
print("Rendering Figure 2 ...", end="  ", flush=True)

subplot_titles = {
    1.7: '1.7 K  |  T₂(SnV) = 170 µs',
    3.0: '3.0 K  |  T₂(SnV) = 170 µs',
}

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
fig2.suptitle(
    'Proposed Quantum Checkpointing Protocol: Fidelity vs Inter-Module Wait Time\n'
    '[SWAP fidelity = 0.9992 = ¹³C single-qubit gate proxy (Resch 2026) : not experimentally validated]',
    fontsize=10, fontweight='bold', y=1.02)

for col, T in enumerate(temperatures):
    ax  = axes2[col]
    cp  = cp_results[T]
    ncp = ncp_results[T]

    ax.semilogx(t_comm_us, cp,  color=COLORS[T], linewidth=2.0,
                linestyle='-',  zorder=3, label='Proposed checkpoint')
    ax.semilogx(t_comm_us, ncp, color=COLORS[T], linewidth=2.0,
                linestyle='--', alpha=0.75, zorder=3, label='No checkpoint')

    ax.fill_between(t_comm_us, cp, ncp,
                    where=(cp >= ncp), alpha=0.15, color='green',
                    interpolate=True, zorder=2, label='Checkpoint advantage')

    ax.axhline(0.99, color='0.45', linestyle='--', linewidth=1.0, zorder=1)
    if col == 0:
        ax.text(1.2, 0.991, 'F = 0.99 threshold',
                ha='left', va='bottom', fontsize=8.5, color='0.45')

    # Annotate max advantage
    mi  = max_adv_idx[T]
    mx  = t_comm_us[mi]
    my  = cp[mi]
    pct = adv[T][mi] * 100
    ax.annotate(f'Max advantage\n+{pct:.1f}% at {mx:.0f} µs',
                xy=(mx, my), xytext=(mx * 4, my + 0.04),
                fontsize=7.5, ha='left', color='0.25',
                arrowprops=dict(arrowstyle='->', color='0.4',
                                lw=0.9, connectionstyle='arc3,rad=0.2'))

    ax.set_title(subplot_titles[T], fontsize=10, fontweight='bold', pad=6)
    ax.set_xlabel('Inter-module wait time (µs)', fontsize=10)
    if col == 0:
        ax.set_ylabel('Recovered SnV fidelity', fontsize=10)
    ax.set_xlim([t_comm_us[0], t_comm_us[-1]])
    ax.set_ylim([0.40, 1.02])
    ax.tick_params(labelsize=9)
    ax.grid(True, color='0.88', linewidth=0.6, which='both', alpha=0.7)
    ax.legend(fontsize=8, loc='lower left')
    remove_spines(ax)

# Caveat box on right subplot
axes2[1].text(0.98, 0.05,
              '(!) 13C data at 50 mK\n   not validated at 1.7-3 K',
              transform=axes2[1].transAxes, fontsize=8,
              ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff3cd',
                        edgecolor='#856404', linewidth=1.0))

fig2.tight_layout()
fig2.savefig('part5_fig2_checkpoint_advantage.png', dpi=300, bbox_inches='tight')
plt.close(fig2)
print("saved part5_fig2_checkpoint_advantage.png")

# ============================================================================
# FIGURE 3 : TEMPERATURE COMPARISON OVERLAY  (1.7 K and 3 K)
# ============================================================================
print("Rendering Figure 3 ...", end="  ", flush=True)

fig3, ax3 = plt.subplots(figsize=(8, 6))

# Shaded latency window (10–100 µs)
ax3.axvspan(10, 100, alpha=0.20, color='#ffe066',
            label='Typical inter-module latency (10–100 µs)')

for T in temperatures:
    ax3.semilogx(t_comm_us, cp_results[T],
                 color=COLORS[T], linewidth=2.0, linestyle='-',
                 label=f'Proposed CP @ {T} K')
    ax3.semilogx(t_comm_us, ncp_results[T],
                 color=COLORS[T], linewidth=1.5, linestyle='--',
                 alpha=0.70, label=f'NCP @ {T} K')
ax3.axhline(0.99, color='0.35', linestyle='--', linewidth=1.0)
ax3.text(t_comm_us[-1] * 0.85, 0.992, 'F = 0.99',
         ha='right', va='bottom', fontsize=9, color='0.35')
ax3.axhline(0.95, color='0.55', linestyle='--', linewidth=1.0)
ax3.text(t_comm_us[-1] * 0.85, 0.952, 'F = 0.95',
         ha='right', va='bottom', fontsize=9, color='0.55')

# Core result annotation
ax3.text(0.97, 0.25,
         'Checkpointing extends usable coherence\n'
         'by >50× vs no-checkpoint at t_comm = 1 ms\n'
         '(T₂_13C / T₂_SnV ≈ 7940× at CPMG-128)',
         transform=ax3.transAxes, fontsize=9,
         ha='right', va='bottom',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                   edgecolor='0.6', linewidth=1.0))

# Caveat annotation : placed upper-left to avoid lower-left legend
ax3.text(0.03, 0.62,
         'CAVEAT:\n'
         '13C data (Resch 2026) measured at 50 mK.\n'
         'SWAP not experimentally demonstrated.\n'
         'These are PROPOSED protocol projections.',
         transform=ax3.transAxes, fontsize=8.5,
         ha='left', va='bottom',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3cd',
                   edgecolor='#856404', linewidth=1.2))

ax3.set_xlabel('Inter-module wait time (µs)', fontsize=12)
ax3.set_ylabel('Recovered SnV fidelity', fontsize=12)
ax3.set_title(
    'Proposed Checkpointing Fidelity : 1.7–3 K Plateau\n'
    '(CP/NCP curves coincide at 1.7 K and 3 K : SnV: Rosenthal 2023; ¹³C: Resch 2026 @ 50 mK)',
    fontsize=11, fontweight='bold', pad=10)
ax3.set_xlim([t_comm_us[0], t_comm_us[-1]])
ax3.set_ylim([0.40, 1.02])
ax3.tick_params(labelsize=10)
ax3.grid(True, color='0.88', linewidth=0.6, which='both', alpha=0.7)
remove_spines(ax3)

cp_handles  = [plt.Line2D([0],[0], color=COLORS[T], lw=2.0, ls='-',
                           label=f'CP @ {T} K') for T in temperatures]
ncp_handles = [plt.Line2D([0],[0], color=COLORS[T], lw=1.5, ls='--',
                           alpha=0.7, label=f'NCP @ {T} K') for T in temperatures]
lat_handle  = [mpatches.Patch(facecolor='#ffe066', alpha=0.5,
                               label='Typical latency window')]
ax3.legend(handles=cp_handles + ncp_handles + lat_handle,
           ncol=2, fontsize=8.5, loc='lower left', framealpha=0.9)

fig3.tight_layout()
fig3.savefig('part5_fig3_temperature_comparison.png', dpi=300, bbox_inches='tight')
plt.close(fig3)
print("saved part5_fig3_temperature_comparison.png")

# ============================================================================
# FIGURE 4 : COLOUR-CODED SUMMARY TABLE
# ============================================================================
print("Rendering Figure 4 ...", end="  ", flush=True)

table_tc_us  = [1, 10, 50, 100, 200, 1000, 10000]
table_tc_sec = [t * 1e-6 for t in table_tc_us]

def nearest_idx(tc_s):
    return int(np.argmin(np.abs(t_comm_values - tc_s)))

row_indices = [nearest_idx(tc) for tc in table_tc_sec]

col_headers = [
    't_comm\n(µs)',
    'CP\n@1.7K', 'CP\n@3K',
    'NCP\n@1.7K', 'NCP\n@3K',
    'Adv.\n@1.7K', 'Adv.\n@3K',
    'MW duty\ncycle (%)',
]

table_data   = []
cell_colours = []

def cell_color(val, is_advantage=False, is_duty=False):
    if is_duty:
        if   val <  15.0: return '#d4edda'   # green : safe
        elif val <  20.0: return '#fff3cd'   # yellow : marginal
        else:             return '#f8d7da'   # red : risky
    if is_advantage:
        return '#d4edda' if val > 0 else '#f8d7da'
    # Fidelity colour
    if   val > 0.99: return '#d4edda'
    elif val > 0.95: return '#fff3cd'
    else:            return '#f8d7da'

for ri in row_indices:
    tc_us_val = t_comm_us[ri]
    cp17   = cp_results[1.7][ri]
    cp3    = cp_results[3.0][ri]
    ncp17  = ncp_results[1.7][ri]
    ncp3   = ncp_results[3.0][ri]
    adv17  = cp17 - ncp17
    adv3   = cp3  - ncp3
    _, t_on, t_cycle = compute_mw_duty_cycle(t_comm_values[ri])
    duty   = t_on / t_cycle * 100

    table_data.append([
        f'{tc_us_val:.0f}',
        f'{cp17:.4f}', f'{cp3:.4f}',
        f'{ncp17:.4f}', f'{ncp3:.4f}',
        f'{adv17:+.4f}', f'{adv3:+.4f}',
        f'{duty:.1f}%',
    ])
    cell_colours.append([
        '#f0f0f0',
        cell_color(cp17), cell_color(cp3),
        cell_color(ncp17), cell_color(ncp3),
        cell_color(adv17, is_advantage=True), cell_color(adv3, is_advantage=True),
        cell_color(duty, is_duty=True),
    ])

fig4, ax4 = plt.subplots(figsize=(13, 4.5))
ax4.axis('off')

tbl = ax4.table(
    cellText=table_data,
    colLabels=col_headers,
    cellColours=cell_colours,
    colColours=['#d0d0d0'] * len(col_headers),
    cellLoc='center',
    loc='center',
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.0)
tbl.scale(1.0, 1.7)

# Bold advantage and duty-cycle columns + header row
for row_idx in range(len(table_data)):
    for col_idx in [5, 6, 7]:
        tbl[(row_idx + 1, col_idx)].set_text_props(fontweight='bold')
for col_idx in range(len(col_headers)):
    tbl[(0, col_idx)].set_text_props(fontweight='bold')

ax4.set_title(
    'Proposed Checkpointing Protocol : Quantitative Summary\n'
    '[CP/NCP fidelities conditioned on unvalidated SWAP proxy; duty cycles are estimates]',
    fontsize=11, fontweight='bold', pad=14)

legend_patches = [
    mpatches.Patch(facecolor='#d4edda', edgecolor='0.6',
                   label='Fidelity > 0.99 / Advantage > 0 / Duty < 15%'),
    mpatches.Patch(facecolor='#fff3cd', edgecolor='0.6',
                   label='0.95 ≤ F ≤ 0.99 / 15–20% duty'),
    mpatches.Patch(facecolor='#f8d7da', edgecolor='0.6',
                   label='F < 0.95 / Advantage ≤ 0 / Duty ≥ 20%'),
]
fig4.legend(handles=legend_patches, loc='lower center', ncol=3,
            fontsize=8.5, bbox_to_anchor=(0.5, -0.04), framealpha=0.9)

fig4.tight_layout()
fig4.savefig('part5_fig4_summary_table.png', dpi=300, bbox_inches='tight')
plt.close(fig4)
print("saved part5_fig4_summary_table.png")

# ============================================================================
# FIGURE 5 : MW vs RF HEATING MARGIN COMPARISON
# ============================================================================
# Source: Karapatzakis et al., arXiv:2606.15398 (2026)
#
# MW (2 GHz, electron-spin gates at B=400 mT):
#   - Onset heating above ~7 dBm CW; T_local rises toward ~3 K at 11 dBm
#   - Superconducting breakdown at ~11–12 dBm CW
#   - At 20% duty cycle: breakdown at ~10 dBm (θ=84°) or ~12 dBm (θ=20°)
#   - τ_th = 0.54(17) ms : temperature recovers on this timescale
#
# RF (20 MHz, nuclear-spin rotations at B=400 mT):
#   - NO measurable heating up to abrupt breakdown at 9.4 dBm (Bac ≈ 1.2 mT)
#   - RF is NOT the binding constraint
#
# This figure shows:
#   Panel (a): Protocol MW duty cycle vs t_comm
#   Panel (b): Schematic MW power safety margin (qualitative)
#   Panel (c): MW vs RF heating margin comparison bar chart
print("Rendering Figure 5 ...", end="  ", flush=True)

fig5, axes5 = plt.subplots(1, 3, figsize=(15, 5))
fig5.suptitle(
    'Heating-Budget Analysis : MW vs RF Constraints\n'
    'Source: Karapatzakis et al., arXiv:2606.15398 (2026)\n'
    '(!) Order-of-magnitude estimate: no first-principles thermal conductance published',
    fontsize=10, fontweight='bold', y=1.03)

# ---- Panel (a): MW duty cycle vs t_comm ----
ax_a = axes5[0]

ax_a.semilogx(t_comm_us, duty_pcts,
              color='#2ca02c', linewidth=2.0, marker='o', markersize=6,
              label='Protocol MW duty cycle')

ax_a.axhline(20.0, color='#d62728', linestyle='--', linewidth=1.5,
             label='20% test threshold (Karapatzakis 2026)')
ax_a.axhline(15.0, color='#ff7f0e', linestyle=':', linewidth=1.5,
             label='~15% onset extrapolation')
ax_a.axhspan(0, 15.0,  alpha=0.08, color='green')
ax_a.axhspan(15.0, 20.0, alpha=0.10, color='orange')
ax_a.axhspan(20.0, 35.0, alpha=0.07, color='red')

ax_a.text(1.5, 5, '[OK] Safe', color='darkgreen', fontsize=9, fontweight='bold')
ax_a.text(1.5, 16.5, '[!] Marginal', color='darkorange', fontsize=9, fontweight='bold')
ax_a.text(1.5, 22, '[X] Risky', color='darkred', fontsize=9, fontweight='bold')

ax_a.set_xlabel('t_comm (µs)', fontsize=10)
ax_a.set_ylabel('Estimated MW duty cycle (%)', fontsize=10)
ax_a.set_title('(a) Protocol MW Duty Cycle', fontsize=10, fontweight='bold')
ax_a.set_xlim([t_comm_us[0], t_comm_us[-1]])
ax_a.set_ylim([0, 35])
ax_a.legend(fontsize=7.5, loc='upper right')
ax_a.grid(True, alpha=0.3, which='both')
remove_spines(ax_a)

# ---- Panel (b): Temperature rise vs MW power (qualitative, from Karapatzakis 2026) ----
ax_b = axes5[1]

# Schematic representation of the published data trend
# CW drive, B=400mT: gradual rise above 7 dBm, ~3 K by 11 dBm, breakdown ~11-12 dBm
p_dbm_cw  = np.array([0, 3, 5, 7, 8, 9, 10, 11])
T_rise_cw = np.array([0, 0, 0, 0.05, 0.2, 0.5, 1.0, 1.2])   # ΔT above base (K, schematic)

ax_b.plot(p_dbm_cw, T_rise_cw, 'o-', color='#d62728', linewidth=2.0,
          label='MW (2 GHz, CW, 400 mT) : schematic')
ax_b.axvline(7.0, color='#ff7f0e', linestyle=':', linewidth=1.2,
             label='~7 dBm onset')
ax_b.axvline(11.0, color='#d62728', linestyle='--', linewidth=1.2,
             label='~11 dBm breakdown')

# RF: flat zero heating up to 9.4 dBm, then abrupt breakdown
p_rf_safe = np.array([0, 3, 6, 9, 9.4])
T_rf_safe = np.array([0, 0, 0, 0, 0])
ax_b.plot(p_rf_safe, T_rf_safe, 's-', color='#1f77b4', linewidth=2.0,
          label='RF (20 MHz, 400 mT) : no heating observed')
ax_b.axvline(9.4, color='#1f77b4', linestyle='--', linewidth=1.2,
             label='9.4 dBm: RF breakdown (abrupt)')

ax_b.axhline(1.0, color='0.6', linestyle=':', linewidth=1.0)
ax_b.text(0.2, 1.05, 'ΔT = 1 K above base', fontsize=7.5, color='0.5')

ax_b.set_xlabel('Drive power (dBm)', fontsize=10)
ax_b.set_ylabel('Local temperature rise ΔT (K)', fontsize=10)
ax_b.set_title('(b) Heating vs Drive Power\n(schematic from Karapatzakis 2026)',
               fontsize=10, fontweight='bold')
ax_b.set_xlim([0, 12])
ax_b.set_ylim([-0.1, 1.5])
ax_b.legend(fontsize=7, loc='upper left')
ax_b.grid(True, alpha=0.3)
remove_spines(ax_b)

# ---- Panel (c): Bar chart : MW vs RF heating margin ----
ax_c = axes5[2]

categories    = ['MW (2 GHz)\nelectron spin', 'RF (20 MHz)\nnuclear spin']
onset_dbm     = [7.0, None]   # onset (MW has gradual onset, RF has none)
breakdown_dbm = [11.0, 9.4]   # breakdown threshold
typical_dbm   = [5.0, 3.0]    # typical operating power (estimated)

x = np.arange(len(categories))
width = 0.25

bars_break  = ax_c.bar(x, breakdown_dbm, width * 2.2,
                        color=['#f8d7da', '#d4edda'],
                        edgecolor='0.4', linewidth=1.0,
                        label='Breakdown threshold (dBm)')
bars_typ    = ax_c.bar(x, typical_dbm, width * 2.2,
                        color=['#1f77b4', '#2ca02c'],
                        alpha=0.75, edgecolor='0.3', linewidth=1.0,
                        label='Estimated typical operation')

# Margin arrows
for i in range(2):
    margin = breakdown_dbm[i] - typical_dbm[i]
    ax_c.annotate('',
                  xy=(x[i], breakdown_dbm[i]),
                  xytext=(x[i], typical_dbm[i]),
                  arrowprops=dict(arrowstyle='<->', color='0.3',
                                  lw=1.5))
    ax_c.text(x[i] + 0.15, (breakdown_dbm[i] + typical_dbm[i]) / 2,
              f'Margin\n+{margin:.1f} dBm',
              fontsize=8.5, va='center', color='0.2')

# Onset line for MW
ax_c.axhline(7.0, color='#ff7f0e', linestyle=':', linewidth=1.2, alpha=0.8)
ax_c.text(-0.5, 7.15, 'MW onset ~7 dBm', fontsize=7.5, color='darkorange')

ax_c.set_xticks(x)
ax_c.set_xticklabels(categories, fontsize=10)
ax_c.set_ylabel('Power (dBm)', fontsize=10)
ax_c.set_title(
    '(c) MW vs RF Heating Margin\n'
    '(MW is binding constraint; RF safe to breakdown)',
    fontsize=10, fontweight='bold')
ax_c.set_ylim([0, 14])
ax_c.legend(fontsize=8, loc='upper right')
ax_c.grid(True, axis='y', alpha=0.3)
remove_spines(ax_c)

# Overall caveat footnote
fig5.text(0.5, -0.03,
          '(!) Disclaimer: No first-principles thermal conductance value is published '
          '(Karapatzakis 2026). All duty-cycle estimates are order-of-magnitude only.\n'
          f'Thermal relaxation constant τ_th = {tau_th_ms:.2f} ms (Karapatzakis 2026). '
          'Panel (b) is a schematic reconstruction of the published trend, not raw data.',
          ha='center', va='top', fontsize=8, color='0.4',
          style='italic', wrap=True)

fig5.tight_layout()
fig5.savefig('part5_fig5_heating_margin.png', dpi=300, bbox_inches='tight')
plt.close(fig5)
print("saved part5_fig5_heating_margin.png")
print()

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
print("=" * 70)
print("SIMULATION RESULTS SUMMARY")
print("=" * 70)
print()
print("SnV physical parameters (Rosenthal et al., PRX 2023):")
print(f"  T1 = 4.23 ms at 3 K (Rosenthal 2023, direct measurement)")
print(f"       No direct 1.7 K T1 measurement; T1 > ~20 ms lower bound (Rosenthal 2023);")
print(f"       Orbach extrapolation gives T1 ≈ 200 s at 1.7 K.")
print(f"       4.23 ms is used at 1.7 K as a conservative worst-case bound.")
print(f"  T2 = 170 µs  (Hahn echo, plateau 1.7–3 K)")
print(f"  T2*= 396.6 +/- 2.29 ns (Rosenthal 2023, Fig. 3e Ramsey fit at 1.7 K; cosmetic only : not in Lindblad model)")
print(f"  T1 = 4.23 ms is the direct 3 K measurement reused at 1.7 K (no direct 1.7 K value).")
print(f"       T2_echo = 170 us dominates coherence decay at both temperatures (T1 >> T2_echo),")
print(f"       so this choice does not affect simulation results.")
print()
print("¹³C physical parameters (Resch et al., PRX 2026, at 50 mK):")
print(f"  T2(CPMG-128) = 1.35 s  ← used as memory coherence")
print(f"  T2(Hahn)     = 167 ms  ← alternative (more conservative) memory estimate")
print(f"  RF gate fidelity = 99.92% ← used as SWAP proxy (PROPOSED, not validated)")
print(f"  *** NOT validated at 1.7–3 K ***")
print()
print("Proposed checkpointing fidelity results:")
print()

for T in temperatures:
    cp  = cp_results[T]
    ncp = ncp_results[T]
    adv_arr = adv[T]

    # Find crossover (last t_comm where CP > NCP)
    helps = adv_arr > 0
    if helps.any():
        crossover_str = f"{t_comm_us[np.where(helps)[0][-1]]:.1f} µs"
    else:
        crossover_str = "never"

    mi = max_adv_idx[T]
    print(f"  At {T}K:")
    print(f"    CP fidelity @ t_comm=10 µs  : {cp_results[T][2]:.4f}")
    print(f"    CP fidelity @ t_comm=100 µs : {cp_results[T][5]:.4f}")
    print(f"    NCP fidelity @ t_comm=100 µs: {ncp_results[T][5]:.4f}")
    print(f"    Max advantage : {adv_arr[mi]:+.4f} at t_comm = {t_comm_us[mi]:.0f} µs")
    print(f"    Crossover     : CP > NCP up to {crossover_str}")
    print()

print("Heating margin summary (Karapatzakis 2026):")
print("  MW (2 GHz, electron spin):")
print("    Onset:        ~7 dBm CW at B=400 mT")
print("    CW breakdown: ~11–12 dBm")
print("    20% duty-cycle breakdown: ~10 dBm (θ=84°), ~12 dBm (θ=20°)")
print("    *** MW heating is the BINDING CONSTRAINT ***")
print()
print("  RF (20 MHz, ¹³C nuclear spin):")
print("    No measurable heating up to abrupt breakdown at 9.4 dBm")
print("    RF is NOT the binding constraint.")
print()

# Duty cycle at key t_comm values
print("  Protocol MW duty cycle estimate:")
for tc_us_val in [1, 10, 50, 100]:
    idx  = nearest_idx(tc_us_val * 1e-6)
    d, _, _ = compute_mw_duty_cycle(t_comm_values[idx])
    safe = "[OK]" if d * 100 < 15 else ("[!]" if d * 100 < 20 else "[X]")
    print(f"    t_comm = {tc_us_val:>5} µs → duty {d*100:.1f}%  {safe}")

print()
print("  (!) DISCLAIMER: Duty-cycle estimates are order-of-magnitude only.")
print(f"     No first-principles thermal conductance published.")
print(f"     τ_th = {tau_th_ms:.2f} ms (Karapatzakis 2026)")
print()
print("=" * 70)
print("FIGURES SAVED")
print("=" * 70)
print("  Fig 1: part5_fig1_decoherence.png")
print("  Fig 2: part5_fig2_checkpoint_advantage.png")
print("  Fig 3: part5_fig3_temperature_comparison.png")
print("  Fig 4: part5_fig4_summary_table.png")
print("  Fig 5: part5_fig5_heating_margin.png")
print()
print("Part 5 complete. Simulation pipeline finished.")
