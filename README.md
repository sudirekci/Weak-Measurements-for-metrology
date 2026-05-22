# Extending the dynamic range in quantum frequency estimation with sequential weak measurements

Code accompanying the paper:

**Su Direkci, Manuel Endres, and Tuvia Gefen**, *Extending the dynamic range in quantum frequency estimation with sequential weak measurements* (2025).
[DOI: 10.1103/2hyn-cqlq](https://doi.org/10.1103/2hyn-cqlq) · [arXiv:2509.01474](https://arxiv.org/abs/2509.01474)

## Overview

This repository contains the numerical simulations used to study sequential weak measurement protocols for Bayesian frequency estimation with coherent spin states. The protocols are designed to extend the dynamic range of Ramsey-type interferometry beyond the standard phase-slip limit (`δω·T > π`), while maintaining asymptotically optimal sensitivity in the large-bandwidth regime.

Two protocols are analyzed:

- **Weak-only protocol** — a sequence of weak `σx` measurements on each sensing qubit during free evolution.
- **Weak-with-strong protocol** — the same sequence, with the final weak measurement replaced by a projective measurement. This protocol asymptotically saturates the quantum Fisher information bound `I_Q = 4NT²` and outperforms cascaded schemes for an arbitrarily large bandwidth.

The code reproduces the main results of the paper, including:

- Classical Fisher information (CFI) for both protocols as a function of interrogation time, measurement strength, and number of qubits.
- Bayesian mean squared error (BMSE) under maximum-likelihood estimation, including the SNR threshold effect.
- Comparison with the optimal classical interferometer (OCI) and the cascaded protocol.
- Analysis with imperfect (noisy) ancilla measurements.

  Examples can be found in sample_fi.m.

## Citation

```bibtex
@article{Direkci2025WeakMeas,
  title   = {Extending the dynamic range in quantum frequency estimation with sequential weak measurements},
  author  = {Direkci, Su and Endres, Manuel and Gefen, Tuvia},
  year    = {2025},
  doi     = {10.1103/2hyn-cqlq},
  eprint  = {2509.01474},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph}
}
```

## Contact

For questions about the code or paper, please contact via email: sdirekci@caltech.edu.
