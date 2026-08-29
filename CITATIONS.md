# CITATIONS.md — References and Parameter Provenance

This repository contains a Lindblad simulation of a **proposed** electron–nuclear spin
checkpointing protocol for tin-vacancy (SnV⁻) diamond qubits at 1.7–3 K, together with an
order-of-magnitude thermal (heating) budget analysis of the MW/RF control lines.

**Every physical parameter used in the code is traceable to the primary sources below.
Every value that is an assumption (not a measurement) is labelled as such.**

---

## 1. Primary data sources (used directly in the code)

### [R1] SnV electron-spin parameters
E. I. Rosenthal, C. P. Anderson, H. C. Kleidermacher, A. J. Stein, H. Lee, J. Grzesik,
G. Scuri, A. E. Rugar, D. Riedel, S. Aghaeimeibodi, G. H. Ahn, K. Van Gasse, J. Vučković,
"Microwave Spin Control of a Tin-Vacancy Qubit in Diamond,"
*Phys. Rev. X* **13**, 031022 (2023). DOI: [10.1103/PhysRevX.13.031022](https://doi.org/10.1103/PhysRevX.13.031022) — arXiv:2306.13199

Values taken: ground-state splitting Δg/2π = 903.0 GHz; T₁ = 4.23 ± 1.37 ms @ 3 K;
T₁ = 5.22 ± 1.54 µs @ 5 K (documented only); T₂,echo = 170 µs (plateau 1.7–3 K);
T₂* = 392.6 ns (Ramsey; documented, not modelled); MW π-pulse t_π = 48.4 ns.

### [R2] ¹³C nuclear-spin parameters
J. Resch, I. Karapatzakis, M. Elshorbagy, M. Schrodin, P. Fuchs, P. Graßhoff, L. Kussi,
C. Sürgers, C. Popov, C. Becher, W. Wernsdorfer, D. Hunger,
"High-Fidelity Control of a ¹³C Nuclear Spin Coupled to a Tin-Vacancy Center in Diamond,"
*Phys. Rev. X* **16**, 011060 (2026). DOI: [10.1103/bmc6-qvwq](https://doi.org/10.1103/bmc6-qvwq) — arXiv:2509.03354

Values taken: T₂(CPMG-128) = 1.35 s (memory coherence input); T₂(Hahn) = 167 ms
(conservative alternative); T₂* = 1.5 ms; initialization fidelity 99.74(3)%;
¹³C single-qubit RF gate fidelity 99.92(1)% (used as the *proposed SWAP proxy*, see caveats).
⚠️ **All values measured at 50 mK; not validated at 1.7–3 K.**

### [R3] Heating thresholds and thermal time constant
I. Karapatzakis, J. Resch, D. Hunger, W. Wernsdorfer,
"Watching a Superconducting Coplanar Waveguide Heat Up with a Single Color Center,"
arXiv:2606.15398 (2026). DOI: [10.48550/arXiv.2606.15398](https://doi.org/10.48550/arXiv.2606.15398)

Values taken: thermal relaxation constant τ_th = 0.54(17) ms; RF (20 MHz, 400 mT):
no measurable heating up to abrupt breakdown at 9.4 dBm (B_ac ≈ 1.2 mT);
MW (GHz, 400 mT): gradual onset ≳ 7 dBm CW, superconducting breakdown ~11–12 dBm CW,
20%-duty-cycle breakdown ~10 dBm (θ = 84°) / ~12 dBm (θ = 20°).

## 2. Checkpointing / network-memory lineage (conceptual basis)

- H. Bernien *et al.*, "Heralded Entanglement Between Solid-State Qubits Separated by
  Three Metres," *Nature* **497**, 86–90 (2013). DOI: 10.1038/nature12016
- A. Reiserer, N. Kalb, M. S. Blok, K. J. M. van Bemmelen, T. H. Taminiau, R. Hanson,
  D. J. Twitchen, M. Markham, "Robust Quantum-Network Memory Using Decoherence-Protected
  Subspaces of Nuclear Spins," *Phys. Rev. X* **6**, 021040 (2016). DOI: 10.1103/PhysRevX.6.021040
- M. K. Bhaskar *et al.*, "Experimental Demonstration of Memory-Enhanced Quantum
  Communication," *Nature* **580**, 60–64 (2020). DOI: 10.1038/s41586-020-2103-5

> **Attribution note:** the decoherence-protected-subspace network-memory paper is
> **Reiserer et al. 2016**, not Cramer et al. 2016 (*Nat. Commun.* **7**, 11526), which is a
> real but separate work (repeated QEC on an encoded logical qubit). Do not conflate the two.

## 3. Supplementary

- I. Karapatzakis *et al.*, "Microwave Control of the Tin-Vacancy Spin Qubit in Diamond
  with a Superconducting Waveguide," *Phys. Rev. X* **14**, 031036 (2024).
  DOI: 10.1103/PhysRevX.14.031036 — waveguide platform underlying [R2]/[R3].

## 4. Parameter provenance table

| Quantity | Value in code | Source | Caveat |
|---|---|---|---|
| SnV Δg/2π | 903.0 GHz | [R1] | — |
| SnV T₁ @ 3 K | 4.23 ms | [R1] | directly measured |
| SnV T₁ @ 1.7 K | 4.23 ms | [R1] (3 K value) | **conservative worst-case bound**; see §5 |
| SnV T₂,echo | 170 µs | [R1] | flat plateau 1.7–3 K |
| SnV T₂* | 392.6 ns | [R1] | documented only, not modelled |
| SnV t_π (MW) | 48.4 ns | [R1] | duty-cycle estimate only |
| ¹³C T₂ (CPMG-128) | 1.35 s | [R2] | **50 mK only** |
| ¹³C RF gate fidelity | 99.92% | [R2] | **SWAP proxy; SWAP not demonstrated** |
| τ_th | 0.54 ms | [R3] | — |
| RF breakdown | 9.4 dBm | [R3] | no heating below it |
| MW onset / breakdown | ~7 / ~11–12 dBm CW | [R3] | 400 mT |
| t_swap | 10 µs | assumption | not measured |
| t_gates | 1 µs | assumption | — |
| pulse filling in SWAP | 10% | assumption | order-of-magnitude |

## 5. Caveats (read before interpreting any output)

1. **T₁ at 1.7 K.** Rosenthal et al. directly measured T₁ at 3 K (4.23 ms) and 5 K (5.22 µs)
   only. At 1.7 K they report no appreciable relaxation over a 20 ms window (a lower bound)
   and an Orbach-model extrapolation of ≈200 s, contingent on the absence of other
   mechanisms. No direct 1.7 K measurement exists. This code conservatively assigns
   T₁(1.7 K) = 4.23 ms as a worst-case lower bound. Results are insensitive to this choice:
   the pure-dephasing rate is constructed so that total coherence decay reproduces the
   measured T₂,echo = 170 µs, and T₂,echo ≪ 2T₁ under either assignment.
2. **¹³C at 50 mK.** All ¹³C values ([R2]) were measured in a dilution refrigerator at
   50 mK. Whether they hold at 1.7–3 K is unverified; this is flagged as future
   experimental work, not assumed.
3. **SWAP proxy.** `swap_fidelity = 0.9992` is the ¹³C *single-qubit* RF gate fidelity from
   [R2]. No SnV↔¹³C SWAP was demonstrated in that experiment (initialization used
   optical+MW pumping). All checkpoint fidelities are optimistic projections conditioned on
   this proxy.
4. **Heating budget is order-of-magnitude.** No first-principles thermal conductance is
   published ([R3]). Duty-cycle estimates are estimates; figure panel 5(b) is a schematic
   reconstruction of the published trend, **not raw data**.

## 6. Code ↔ source map

| File | Uses |
|---|---|
| `part1_system_setup.py` | none (Hilbert-space setup only) |
| `part2_decoherence.py` | [R1], [R2] |
| `part3_swap.py` | [R2] (proxy) |
| `part4_checkpoint.py` | [R1], [R2], [R3] |
| `part5_figures.py` | [R1], [R2], [R3] |

## 7. BibTeX

```bibtex
@article{Rosenthal2023,
  author  = {Rosenthal, Eric I. and Anderson, Christopher P. and Kleidermacher, Hannah C.
             and Stein, Abigail J. and Lee, Hope and Grzesik, Jakob and Scuri, Giovanni
             and Rugar, Ammon E. and Riedel, David and Aghaeimeibodi, Sattar
             and Ahn, Geun Ho and Van Gasse, Kasper and Vu\v{c}kovi\'c, Jelena},
  title   = {Microwave Spin Control of a Tin-Vacancy Qubit in Diamond},
  journal = {Phys. Rev. X}, volume = {13}, pages = {031022}, year = {2023},
  doi     = {10.1103/PhysRevX.13.031022}, note = {arXiv:2306.13199}}

@article{Resch2026,
  author  = {Resch, Jeremias and Karapatzakis, Ioannis and Elshorbagy, Mohamed
             and Schrodin, Marcel and Fuchs, Philipp and Gra\ss{}hoff, Philipp
             and Kussi, Luis and S\"urgers, Christoph and Popov, Cyril
             and Becher, Christoph and Wernsdorfer, Wolfgang and Hunger, David},
  title   = {High-Fidelity Control of a $^{13}$C Nuclear Spin Coupled to a
             Tin-Vacancy Center in Diamond},
  journal = {Phys. Rev. X}, volume = {16}, pages = {011060}, year = {2026},
  doi     = {10.1103/bmc6-qvwq}, note = {arXiv:2509.03354}}

@misc{Karapatzakis2026,
  author  = {Karapatzakis, Ioannis and Resch, Jeremias and Hunger, David
             and Wernsdorfer, Wolfgang},
  title   = {Watching a Superconducting Coplanar Waveguide Heat Up
             with a Single Color Center},
  year    = {2026}, eprint = {2606.15398},
  archivePrefix = {arXiv}, primaryClass = {quant-ph},
  doi     = {10.48550/arXiv.2606.15398}}

@article{Bernien2013,
  author  = {Bernien, H. and Hensen, B. and Pfaff, W. and Koolstra, G. and Blok, M. S.
             and Robledo, L. and Taminiau, T. H. and Markham, M. and Twitchen, D. J.
             and Morton, J. J. L. and Hanson, R.},
  title   = {Heralded Entanglement Between Solid-State Qubits Separated by Three Metres},
  journal = {Nature}, volume = {497}, pages = {86--90}, year = {2013},
  doi     = {10.1038/nature12016}}

@article{Reiserer2016,
  author  = {Reiserer, Andreas and Kalb, Norbert and Blok, Machiel S.
             and van Bemmelen, Koen J. M. and Taminiau, Tim H. and Hanson, Ronald
             and Twitchen, David J. and Markham, Mark},
  title   = {Robust Quantum-Network Memory Using Decoherence-Protected
             Subspaces of Nuclear Spins},
  journal = {Phys. Rev. X}, volume = {6}, pages = {021040}, year = {2016},
  doi     = {10.1103/PhysRevX.6.021040}, note = {arXiv:1603.01602}}

@article{Bhaskar2020,
  author  = {Bhaskar, Mihir K. and Riedinger, Ralf and Machielse, Bartholomeus
             and Levonian, David S. and Nguyen, Charles T. and Knall, Erik N.
             and Park, Hongkun and Englund, Dirk and Lon\v{c}ar, Marko
             and Sukachev, Dmitry D. and Lukin, Mikhail D.},
  title   = {Experimental Demonstration of Memory-Enhanced Quantum Communication},
  journal = {Nature}, volume = {580}, pages = {60--64}, year = {2020},
  doi     = {10.1038/s41586-020-2103-5}, note = {arXiv:1909.01323}}

@article{Karapatzakis2024,
  author  = {Karapatzakis, Ioannis and Resch, Jeremias and Schrodin, Marcel
             and Fuchs, Philipp and Kieschnick, Michael and Heupel, Julia
             and Kussi, Luis and S\"urgers, Christoph and Popov, Cyril
             and Meijer, Jan and Becher, Christoph and Wernsdorfer, Wolfgang
             and Hunger, David},
  title   = {Microwave Control of the Tin-Vacancy Spin Qubit in Diamond
             with a Superconducting Waveguide},
  journal = {Phys. Rev. X}, volume = {14}, pages = {031036}, year = {2024},
  doi     = {10.1103/PhysRevX.14.031036}}
```
