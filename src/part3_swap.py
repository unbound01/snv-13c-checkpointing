"""
Part 3: SWAP Gate Implementation and Verification
==================================================
Diamond quantum computing architecture — quantum checkpointing protocol.

Qubit 1: SnV electron spin  (compute qubit)
Qubit 2: 13C nuclear spin   (memory qubit)

This file covers:
  - SWAP unitary construction and property verification
  - Noisy SWAP via depolarising channel
  - Forward SWAP:  SnV |+⟩ → 13C,  SnV lands in |0⟩
  - Reverse SWAP:  13C |+⟩ → SnV,  SnV recovers |+⟩ with F≈0.9984
  - Bloch sphere visualisation at each stage

IMPORTANT LABELLING NOTE — READ BEFORE INTERPRETING RESULTS
-------------------------------------------------------------
The "swap_fidelity" parameter (0.9992) is sourced from the single-qubit RF gate
fidelity for ¹³C nuclear spin rotations reported in:
  Resch et al., Phys. Rev. X 16, 011060 (2026)  [arXiv:2509.03354]

This is NOT an experimentally measured or demonstrated SWAP fidelity.
In Resch et al. the ¹³C qubit was initialized via simultaneous optical + MW
pumping — no two-qubit SWAP operation between SnV electron spin and ¹³C was
demonstrated or reported.  The 0.9992 value is used here as a proxy for a
*proposed/unvalidated* SWAP gate, representing an optimistic upper bound on
what a future hyperfine-mediated SWAP might achieve.

Any figure or table derived from this file should be captioned:
  "Proposed checkpointing protocol assuming a SWAP gate fidelity equal to the
   demonstrated ¹³C single-qubit gate fidelity (Resch et al. 2026, 50 mK).
   The SWAP operation itself has not been experimentally demonstrated."

No wait time. No decoherence during SWAP. No time evolution.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, qeye, ket2dm, ptrace, expect,
    sigmax, sigmay, sigmaz, Qobj, Bloch
)

# ============================================================================
# SYSTEM REBUILD  (standalone — consistent with Parts 1 & 2)
# ============================================================================
I2 = qeye(2)

snv_0 = basis(2, 0)
snv_1 = basis(2, 1)
c13_0 = basis(2, 0)
c13_1 = basis(2, 1)

psi_snv_plus = (snv_0 + snv_1).unit()
psi_snv_zero = snv_0
psi_13c_zero = c13_0
psi_13c_plus = (c13_0 + c13_1).unit()

psi_full_init = tensor(psi_snv_plus, psi_13c_zero)
rho_init      = ket2dm(psi_full_init)

# ============================================================================
# PHYSICAL PARAMETERS
# ============================================================================
# swap_fidelity = 0.9992
# SOURCE: Resch et al., PRX 16, 011060 (2026) — ¹³C RF single-qubit gate fidelity
# *** THIS IS NOT AN EXPERIMENTALLY CONFIRMED SWAP GATE FIDELITY ***
# It is used as a proxy for a proposed/unvalidated SnV↔¹³C SWAP operation.
# The 13C data were taken at 50 mK. No SWAP was performed in that experiment.
swap_fidelity = 0.9992   # proposed nuclear-spin single-qubit gate fidelity (Resch 2026)
swap_time     = 10e-6    # 10 µs gate duration (assumed, not measured for a SWAP)

# ============================================================================
# HELPER: PURE-STATE FIDELITY
# ============================================================================
def pure_fidelity(psi_target, rho):
    """Return F = <psi_target|rho|psi_target> as a real float in [0,1]."""
    return float((psi_target.dag() * rho * psi_target).real)

# ============================================================================
# 1. SWAP UNITARY
# ============================================================================
swap_matrix = np.array([[1, 0, 0, 0],
                         [0, 0, 1, 0],
                         [0, 1, 0, 0],
                         [0, 0, 0, 1]], dtype=complex)

U_swap = Qobj(swap_matrix, dims=[[2, 2], [2, 2]])

# ============================================================================
# 2. SWAP PROPERTY VERIFICATION
# ============================================================================
I4 = tensor(I2, I2)

is_unitary    = (U_swap.dag() * U_swap - I4).norm() < 1e-10
is_hermitian  = (U_swap.dag() - U_swap).norm()       < 1e-10
is_involution = (U_swap * U_swap - I4).norm()         < 1e-10

print("=" * 70)
print("PART 3: SWAP GATE IMPLEMENTATION AND VERIFICATION")
print("Diamond Quantum Computing — Checkpointing Protocol")
print("=" * 70)
print()
print("NOTE: 'SWAP fidelity' = 0.9992 is the ¹³C single-qubit RF gate fidelity")
print("  from Resch et al. 2026. No SWAP between SnV and ¹³C was demonstrated.")
print("  All fidelity results below represent a PROPOSED/UNVALIDATED operation.")
print()
print("─── SWAP Matrix Properties ───────────────────────────────")
print(f"  Is SWAP unitary?       (U†U = I)   →  {is_unitary}")
print(f"  Is SWAP hermitian?     (U† = U)    →  {is_hermitian}")
print(f"  Is SWAP self-inverse?  (U² = I)    →  {is_involution}")
print()
assert is_unitary    and is_hermitian and is_involution
print("  All SWAP property assertions passed.")
print()

# ============================================================================
# 3. NOISY SWAP FUNCTION
# ============================================================================
# Depolarising channel:
#   ρ_noisy = F_swap × ρ_ideal + (1 − F_swap) × I/4
#
# Using swap_fidelity = 0.9992 (¹³C single-qubit gate fidelity, Resch 2026)
# as proxy for proposed SWAP — NOT experimentally validated.

rho_mixed = I4 / 4.0

def noisy_swap(rho):
    """
    Apply one noisy SWAP gate (depolarising channel).
    F used = 0.9992 = ¹³C RF gate fidelity (Resch 2026 @ 50 mK).
    Proxy for a proposed SnV↔¹³C SWAP — not experimentally demonstrated.
    """
    rho_ideal = U_swap * rho * U_swap.dag()
    return swap_fidelity * rho_ideal + (1.0 - swap_fidelity) * rho_mixed

# ============================================================================
# 4. FORWARD SWAP:  SnV |+⟩  →  13C
# ============================================================================
rho_after_swap_in = noisy_swap(rho_init)

rho_snv_after_in = ptrace(rho_after_swap_in, 0)
rho_13c_after_in = ptrace(rho_after_swap_in, 1)

fid_snv_after_in = pure_fidelity(psi_snv_zero, rho_snv_after_in)
fid_13c_after_in = pure_fidelity(psi_13c_plus, rho_13c_after_in)

print("─── Forward SWAP (proposed): SnV |+⟩ → ¹³C ──────────────")
print()
print("  SnV reduced density matrix (should be ≈ |0⟩⟨0|):")
print(rho_snv_after_in)
print()
print("  13C reduced density matrix (should be ≈ |+⟩⟨+|):")
print(rho_13c_after_in)
print()
print(f"  SnV fidelity with |0⟩  :  {fid_snv_after_in:.6f}  (expected ≈ {swap_fidelity:.4f})")
print(f"  13C fidelity with |+⟩  :  {fid_13c_after_in:.6f}  (expected ≈ {swap_fidelity:.4f})")
print()

# ============================================================================
# 5. REVERSE SWAP:  13C |+⟩  →  SnV
# ============================================================================
rho_after_swap_out = noisy_swap(rho_after_swap_in)

rho_snv_after_out = ptrace(rho_after_swap_out, 0)
rho_13c_after_out = ptrace(rho_after_swap_out, 1)

fid_snv_after_out = pure_fidelity(psi_snv_plus, rho_snv_after_out)
fid_13c_after_out = pure_fidelity(psi_13c_zero, rho_13c_after_out)

fid_roundtrip_analytical = swap_fidelity ** 2

print("─── Reverse SWAP (proposed): ¹³C |+⟩ → SnV ──────────────")
print()
print("  SnV reduced density matrix (should be ≈ |+⟩⟨+|):")
print(rho_snv_after_out)
print()
print("  13C reduced density matrix (should be ≈ |0⟩⟨0|):")
print(rho_13c_after_out)
print()
print(f"  SnV fidelity with |+⟩  :  {fid_snv_after_out:.6f}  (expected ≈ {fid_roundtrip_analytical:.6f})")
print(f"  13C fidelity with |0⟩  :  {fid_13c_after_out:.6f}  (expected ≈ {fid_roundtrip_analytical:.6f})")
print()

# ============================================================================
# 6. SUMMARY TABLE
# ============================================================================
print("=" * 70)
print("SWAP SUMMARY TABLE  [PROPOSED OPERATION — see caveat above]")
print("=" * 70)
header = f"{'Stage':<18} {'SnV state':<12} {'13C state':<12} {'F_SnV':>8}  {'F_13C':>8}"
print(header)
print("-" * 70)
print(f"{'Initial':<18} {'|+⟩':<12} {'|0⟩':<12} {1.0:>8.4f}  {1.0:>8.4f}")
print(f"{'After SWAP in':<18} {'|0⟩':<12} {'|+⟩':<12} {fid_snv_after_in:>8.4f}  {fid_13c_after_in:>8.4f}")
print(f"{'After SWAP out':<18} {'|+⟩':<12} {'|0⟩':<12} {fid_snv_after_out:>8.4f}  {fid_13c_after_out:>8.4f}")
print()

# ============================================================================
# 7. INTERPRETATION
# ============================================================================
fid_roundtrip_exact = (swap_fidelity**2 + (1 - swap_fidelity**2) * 0.25)

print("─── Interpretation ───────────────────────────────────────────────────")
print(f"  Forward SWAP fidelity          :  {fid_snv_after_in:.6f}")
print(f"  Reverse SWAP fidelity          :  {fid_snv_after_out:.6f}")
print(f"  Round-trip fidelity            :  {fid_snv_after_out:.6f}")
print(f"  Analytical prediction  F²      :  {fid_roundtrip_analytical:.6f}")
print(f"  Exact depolarising model pred  :  {fid_roundtrip_exact:.6f}")
print(f"    (= F² + (1−F²)×0.25, includes mixed-state floor)")
print(f"  Simulation vs exact prediction :  {abs(fid_snv_after_out - fid_roundtrip_exact):.2e}")
print()
print("  REMINDER: F = 0.9992 = ¹³C single-qubit gate fidelity (Resch 2026).")
print("  A SnV↔¹³C SWAP was NOT demonstrated in that experiment.")
print("  Round-trip cost of 0.16% is an OPTIMISTIC LOWER BOUND.")
print()

# ============================================================================
# 8. BLOCH SPHERE PLOTS  (3 rows × 2 columns)
# ============================================================================
# Row 1: Initial state        — SnV=|+⟩,  13C=|0⟩
# Row 2: After forward SWAP   — SnV≈|0⟩,  13C≈|+⟩
# Row 3: After reverse SWAP   — SnV≈|+⟩,  13C≈|0⟩
#
# FIX: Use fig.add_subplot(..., projection='3d') individually (not plt.subplots
# with subplot_kw), then pass the axes to Bloch(axes=ax).  The plt.subplots
# approach creates axes objects that lack the 3D projection attributes that
# QuTiP's Bloch.make_sphere() queries, producing blank output.

snv_color = '#1f77b4'
c13_color = '#d62728'

snv_rhos = [ptrace(rho_init, 0),
            rho_snv_after_in,
            rho_snv_after_out]

c13_rhos = [ptrace(rho_init, 1),
            rho_13c_after_in,
            rho_13c_after_out]

snv_titles = [
    f"SnV : |+⟩  (F={1.0:.4f})",
    f"SnV : ≈|0⟩  (F={fid_snv_after_in:.4f})",
    f"SnV : ≈|+⟩  (F={fid_snv_after_out:.4f})",
]
c13_titles = [
    f"¹³C : |0⟩  (F={1.0:.4f})",
    f"¹³C : ≈|+⟩  (F={fid_13c_after_in:.4f})",
    f"¹³C : ≈|0⟩  (F={fid_13c_after_out:.4f})",
]

fig = plt.figure(figsize=(10, 14))
fig.suptitle(
    "Part 3 : Proposed SWAP Gate: Qubit States at Each Stage\n"
    "[F=0.9992 is ¹³C single-qubit gate fidelity : SWAP not demonstrated experimentally]",
    fontsize=12, fontweight='bold', y=0.99)

for row in range(3):
    # SnV (left column) — created via add_subplot, not plt.subplots
    ax_snv = fig.add_subplot(3, 2, row * 2 + 1, projection='3d')
    b_snv = Bloch(fig=fig, axes=ax_snv)
    b_snv.add_states(snv_rhos[row])
    b_snv.vector_color = [snv_color]
    b_snv.point_color  = [snv_color]
    b_snv.frame_alpha  = 0.1
    b_snv.render()
    ax_snv.set_title(snv_titles[row], fontsize=9, pad=10)

    # 13C (right column)
    ax_13c = fig.add_subplot(3, 2, row * 2 + 2, projection='3d')
    b_13c = Bloch(fig=fig, axes=ax_13c)
    b_13c.add_states(c13_rhos[row])
    b_13c.vector_color = [c13_color]
    b_13c.point_color  = [c13_color]
    b_13c.frame_alpha  = 0.1
    b_13c.render()
    ax_13c.set_title(c13_titles[row], fontsize=9, pad=10)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('part3_swap_bloch.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print("Plot saved: part3_swap_bloch.png")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("=" * 70)
print("PHYSICAL SUMMARY  [PROPOSED OPERATION CAVEAT APPLIES]")
print("=" * 70)
print()
print(f"  Nuclear-spin single-qubit gate fidelity used as SWAP proxy : {swap_fidelity}")
print(f"  Source: Resch et al., PRX 16, 011060 (2026)  [at 50 mK]")
print(f"  SWAP gate duration (assumed)                              : {swap_time*1e6:.0f} µs")
print()
print("  Depolarising noise model:")
print("    ρ_noisy = F_swap × ρ_ideal + (1−F_swap) × I/4")
print()
print("  Forward SWAP transfers SnV |+⟩ into ¹³C with F ≈ 0.9992")
print("  Reverse SWAP recovers SnV |+⟩ from ¹³C  with F ≈ 0.9984")
print("  Round-trip loss = 1 − 0.9992² = 0.0016  (0.16%)")
print()
print("  *** IMPORTANT: This is an optimistic bound. A real SWAP would")
print("      likely incur additional overhead from the two-qubit interaction,")
print("      motional decoherence, and cross-talk. ***")
print()
print("Part 3 complete. Ready for Part 4: Full Checkpointing Protocol.")
