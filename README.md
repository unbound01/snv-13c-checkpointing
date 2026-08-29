# SnV–¹³C Checkpointing Simulation

A Lindblad-equation simulation of an electron–nuclear spin checkpointing protocol for tin-vacancy (SnV) diamond qubits at 1.7–3K, including a heating-budget feasibility analysis grounded in published experimental data.

## What this is

This repository models a proposed two-tier quantum memory scheme for SnV diamond qubits: a fast electron-spin "compute" qubit and a long-lived ¹³C nuclear-spin "memory" qubit, connected by a hyperfine SWAP operation. The idea is to protect the electron spin's quantum state during long communication delays by temporarily storing it in the more stable nuclear spin.

The simulation quantifies:
- Decoherence of both qubits at 1.7K and 3K, using directly measured SnV parameters
- Fidelity of a proposed checkpointing protocol (SWAP-in → wait → SWAP-out) versus no checkpointing
- Whether the microwave (electron-spin) and RF (nuclear-spin) control pulses required for this protocol stay within measured thermal safety limits

## What this is NOT

This is a **feasibility study built on real experimental parameters**, not a demonstrated device or a claim that the full protocol has been built. Specifically:

- The ¹³C nuclear spin coherence times used here (T2 = 1.35s via CPMG-128) were measured at **50mK in a dilution refrigerator**, not at the 1.7–3K operating point used for the SnV electron spin. Whether this coherence survives at 1.7–3K is an open question this simulation does not answer.
- **No SWAP gate between SnV and ¹³C has been experimentally demonstrated anywhere in the cited literature.** The SWAP fidelity used here (0.9992) is the ¹³C single-qubit RF gate fidelity from Resch et al., used as a proxy upper bound for a proposed, not-yet-built SWAP operation.
- The heating-budget duty-cycle analysis is an order-of-magnitude estimate. No published first-principles thermal conductance value exists for this system; the estimate interpolates from measured power-vs-heating data points.

These caveats are stated explicitly in the relevant figures and are not incidental — they define the actual scope of what this work claims.

## Physical parameters and sources

| Quantity | Value | Source |
|---|---|---|
| SnV ground-state splitting | Δg/2π = 903.0 ± 0.7 GHz | Rosenthal et al., *PRX* 13, 031022 (2023) |
| SnV T1 | 4.23 ± 1.37 ms (measured at 3K) | Rosenthal et al. 2023 |
| SnV T2 (Hahn echo) | 170.0 ± 2.8 µs (plateau, 1.7–3K) | Rosenthal et al. 2023 |
| SnV T2* (Ramsey) | 396.6 ± 2.29 ns (at 1.7K) | Rosenthal et al. 2023 |
| ¹³C T2 (128-pulse CPMG) | 1.35 ± 0.03 s (at 50mK) | Resch et al., *PRX* 16, 011060 (2026) |
| ¹³C initialization fidelity | 99.74 ± 0.03% | Resch et al. 2026 |
| ¹³C RF gate fidelity | 99.92 ± 0.01% | Resch et al. 2026 |
| MW/RF heating thresholds | see Fig. 5 | Karapatzakis et al., arXiv:2606.15398 (2026) |

Full reference list with DOIs and arXiv IDs in [`CITATIONS.md`](CITATIONS.md).

Note: SnV T1 has no direct measurement at 1.7K; the 3K value is conservatively reused (see `part2_decoherence.py` and Fig. 1 caption). This does not materially affect results since T2(Hahn echo) is the dominant coherence-limiting timescale at both temperatures.

## Repository structure

```
src/
  part1_system_setup.py      Hilbert space, state initialization, Bloch sphere verification
  part2_decoherence.py       Lindblad decoherence simulation at 1.7K and 3K
  part3_swap.py              Noisy SWAP gate model and verification
  part4_checkpoint.py        Full checkpointing protocol simulation
  part5_figures.py           Publication figure generation, including heating-budget analysis
figures/
  All generated PNG outputs
CITATIONS.md                 Full reference list with DOIs/arXiv IDs
requirements.txt             Python dependencies
LICENSE                      MIT
```

## Running the simulation

```bash
pip install -r requirements.txt
python src/part1_system_setup.py
python src/part2_decoherence.py
python src/part3_swap.py
python src/part4_checkpoint.py
python src/part5_figures.py
```

Scripts are standalone and can be run independently; each rebuilds its own system state. Figures are written to the working directory / `figures/`.

## Key results

- Proposed checkpointing maintains fidelity ≈ 0.968 across inter-module wait times from 1µs to 10ms; the no-checkpoint baseline collapses to ≈0.50 by ~1ms.
- Estimated microwave duty cycle for the protocol stays below the measured 20% breakdown-test threshold across all simulated communication times, with comfortable margin (≤13.3% even at the shortest, most demanding case).
- RF (nuclear-spin) driving shows no measurable heating in the source data up to its own breakdown threshold — the electron-spin MW step, not the RF step, is the binding thermal constraint.

All numeric results are conditioned on the SWAP-proxy and 50mK caveats above.

## Citation

If you use this code, please cite the corresponding paper (link/DOI to be added on publication) and the original experimental sources listed in `CITATIONS.md`.

## License

MIT — see [`LICENSE`](LICENSE).
