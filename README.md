# Braid Monodromy Rigidity on One-Dimensional Symmetric Cyclic Prym Slices: Adjoint $\mathfrak{sl}_4$-Dichotomy, Ramanujan-Sum Trace Averages, and Length-Three Cyclotomic Separation

[![verify](https://github.com/Li-Hongmin/cyclic-prym-braid-monodromy/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/cyclic-prym-braid-monodromy/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054607-blue)](https://doi.org/10.5281/zenodo.24054607)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: on one-dimensional symmetric slices of cyclic Prym families, braid monodromy is rigid, with an adjoint sl₄ dichotomy, Ramanujan-sum trace averages and length-three cyclotomic separation.

- Manuscript on rStation: <https://rstation.ai/problems/hodge_conjecture/papers/114/v4> (version 4)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054607](https://doi.org/10.5281/zenodo.24054607) (all versions); this repository mirrors version 4, DOI [10.5281/zenodo.24054633](https://doi.org/10.5281/zenodo.24054633)

## Contents

- `paper/paper.md`, `paper/cyclic-prym-braid-monodromy.pdf` — the manuscript (version 4).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Slice monodromy computations | `python verify_slice_monodromy.py` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
sha256sum -c CONTENTS.sha256      # or: shasum -a 256 -c CONTENTS.sha256
pip install sympy==1.14.0
python verify_slice_monodromy.py
```

Requirements: Python 3.9+. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026cyclicprymbraidmonodromy,
  author    = {Li, Hongmin},
  title     = {Braid Monodromy Rigidity on One-Dimensional Symmetric Cyclic Prym Slices: Adjoint $\mathfrak{sl}_4$-Dichotomy, Ramanujan-Sum Trace Averages, and Length-Three Cyclotomic Separation},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054607},
  url       = {https://doi.org/10.5281/zenodo.24054607}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
