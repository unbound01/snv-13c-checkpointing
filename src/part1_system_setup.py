"""
Part 1: Two-Qubit System Setup and State Verification
======================================================
Diamond quantum computing architecture : quantum checkpointing protocol.

Qubit 1: SnV electron spin  (compute qubit)
Qubit 2: 13C nuclear spin   (memory qubit)

This file ONLY covers:
  - Hilbert space definition
  - Operator construction on the full tensor-product space
  - State initialisation
  - State verification (density matrices, fidelities, Bloch vectors)
  - Bloch sphere plots

No decoherence. No SWAP. No time evolution.

BLOCH SPHERE RENDERING FIX
-----------------------------
The previous version used plt.subplots(subplot_kw=dict(projection='3d'))
which creates a Figure with axes objects that do not carry the full 3D
projection state expected by QuTiP's Bloch.make_sphere().  This caused
blank white output.  The fix is to create the figure with plt.figure()
and add each 3D subplot individually via fig.add_subplot(..., projection='3d'),
then pass those axes to Bloch(axes=ax).  This is the pattern that works
reliably across QuTiP 4.x and 5.x.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, qeye, sigmax, sigmay, sigmaz,
    ket2dm, ptrace, expect, Bloch
)

# ============================================================================
# 1. HILBERT SPACE DEFINITION
# ============================================================================
I2 = qeye(2)

snv_0 = basis(2, 0)
snv_1 = basis(2, 1)
c13_0 = basis(2, 0)
c13_1 = basis(2, 1)

print("=" * 60)
print("PART 1: TWO-QUBIT SYSTEM SETUP")
print("Diamond Quantum Computing : Checkpointing Protocol")
print("=" * 60)
print()
print("Hilbert space: H_SnV (dim=2) ⊗ H_13C (dim=2)  →  dim=4")
print()

# ============================================================================
# 2. OPERATORS ON THE FULL HILBERT SPACE
# ============================================================================
Sx_snv = tensor(sigmax(), I2)
Sy_snv = tensor(sigmay(), I2)
Sz_snv = tensor(sigmaz(), I2)

Sx_13c = tensor(I2, sigmax())
Sy_13c = tensor(I2, sigmay())
Sz_13c = tensor(I2, sigmaz())

print("Operators defined on full Hilbert space (4×4 matrices):")
print("  SnV: Sx_snv, Sy_snv, Sz_snv  =  σ_i ⊗ I")
print("  13C: Sx_13c, Sy_13c, Sz_13c  =  I ⊗ σ_i")
print()

# ============================================================================
# 3. STATE INITIALISATION
# ============================================================================
psi_snv  = (snv_0 + snv_1).unit()   # |+>_SnV
psi_13c  = c13_0                     # |0>_13C
psi_full = tensor(psi_snv, psi_13c)
rho_full = ket2dm(psi_full)

print("Initial states:")
print("  SnV  : |+> = (|0> + |1>) / sqrt(2)")
print("  13C  : |0>")
print("  Full : |+>_SnV ⊗ |0>_13C")
print()

# ============================================================================
# 4. VERIFICATION
# ============================================================================
rho_snv = ptrace(rho_full, 0)
rho_13c = ptrace(rho_full, 1)

print("-" * 60)
print("REDUCED DENSITY MATRICES")
print("-" * 60)
print()
print("SnV reduced density matrix  ρ_SnV:")
print(rho_snv)
print()
print("13C reduced density matrix  ρ_13C:")
print(rho_13c)
print()

fid_snv = float((psi_snv.dag() * rho_snv * psi_snv).real)
fid_13c = float((psi_13c.dag() * rho_13c * psi_13c).real)

print("-" * 60)
print("FIDELITY CHECKS")
print("-" * 60)
print()
print(f"  SnV fidelity with |+>  :  {fid_snv:.6f}  (expected 1.000000)")
print(f"  13C fidelity with |0>  :  {fid_13c:.6f}  (expected 1.000000)")
print()

rx_snv = expect(sigmax(), rho_snv)
ry_snv = expect(sigmay(), rho_snv)
rz_snv = expect(sigmaz(), rho_snv)
rx_13c = expect(sigmax(), rho_13c)
ry_13c = expect(sigmay(), rho_13c)
rz_13c = expect(sigmaz(), rho_13c)

print("-" * 60)
print("BLOCH VECTORS")
print("-" * 60)
print()
print(f"  SnV  [{rx_snv:+.6f},  {ry_snv:+.6f},  {rz_snv:+.6f}]  (expected [+1, 0, 0])")
print(f"  13C  [{rx_13c:+.6f},  {ry_13c:+.6f},  {rz_13c:+.6f}]  (expected [0, 0, +1])")
print()

assert abs(fid_snv - 1.0) < 1e-9, f"SnV fidelity is not 1! Got {fid_snv}"
assert abs(fid_13c - 1.0) < 1e-9, f"13C fidelity is not 1! Got {fid_13c}"
assert abs(rx_snv - 1.0)  < 1e-9, f"SnV Bloch x should be 1! Got {rx_snv}"
assert abs(rz_13c - 1.0)  < 1e-9, f"13C Bloch z should be 1! Got {rz_13c}"
print("  All assertions passed : states correctly initialised.")
print()

# ============================================================================
# 5. BLOCH SPHERE PLOTS
# ============================================================================
# FIX: Create figure with plt.figure(), then add each 3D subplot individually
# via fig.add_subplot(..., projection='3d').  Do NOT use plt.subplots with
# subplot_kw=dict(projection='3d') : that approach produces blank white output
# when the axes are passed to QuTiP's Bloch(axes=ax) because the axes objects
# created by plt.subplots do not retain the required 3D projection attributes
# in all QuTiP versions.

fig = plt.figure(figsize=(10, 5))
fig.suptitle("Part 1 : Initial Qubit States on Bloch Spheres",
             fontsize=14, fontweight='bold', y=1.01)

# SnV : left subplot
ax_snv = fig.add_subplot(1, 2, 1, projection='3d')
b_snv  = Bloch(fig=fig, axes=ax_snv)
b_snv.add_states(rho_snv)
b_snv.vector_color = ['#1f77b4']
b_snv.point_color  = ['#1f77b4']
b_snv.frame_alpha  = 0.1
b_snv.render()
ax_snv.set_title("SnV electron spin\n|+⟩ = (|0⟩+|1⟩)/√2",
                 fontsize=11, pad=12)

# 13C : right subplot
ax_13c = fig.add_subplot(1, 2, 2, projection='3d')
b_13c  = Bloch(fig=fig, axes=ax_13c)
b_13c.add_states(rho_13c)
b_13c.vector_color = ['#d62728']
b_13c.point_color  = ['#d62728']
b_13c.frame_alpha  = 0.1
b_13c.render()
ax_13c.set_title("¹³C nuclear spin\n|0⟩ (north pole)",
                 fontsize=11, pad=12)

plt.tight_layout()
plt.savefig("part1_bloch_spheres.png", dpi=150, bbox_inches='tight')
plt.close(fig)
print("-" * 60)
print("Bloch sphere plot saved as:  part1_bloch_spheres.png")
print("-" * 60)
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Full system Hilbert space dimension : {rho_full.shape[0]}")
print(f"  SnV state                           : |+>")
print(f"  13C state                           : |0>")
print(f"  SnV fidelity with |+>               : {fid_snv:.6f}")
print(f"  13C fidelity with |0>               : {fid_13c:.6f}")
print(f"  SnV Bloch vector                    : [{rx_snv:+.4f}, {ry_snv:+.4f}, {rz_snv:+.4f}]")
print(f"  13C Bloch vector                    : [{rx_13c:+.4f}, {ry_13c:+.4f}, {rz_13c:+.4f}]")
print(f"  Entanglement (product state)        : None : separable by construction")
print()
print("Part 1 complete. Ready for Part 2: Decoherence.")
