# IR-limits

This project tests the various Infrared (IR) limits necessary for the construction of an N3LO subtraction scheme.

Note, that the repository relies on dependencies like `Stripper` and `Otter` which are not publicly available.

---

## Features

The repository has 3 subdirectories

- `form`: contains programs for analytic manipulations. Specifically it is designed to
  * square soft currents
  * perform UV- and IR-renormalization
- `latex`: contains the LaTeX source files for the analytic results of the IR limits and readable form and explains the used conventions.
- `numerics`: This contains programs to test and validate the analytic expressions numerically.

---

## Installation

Clone the repository

```bash
git clone https://github.com/tomtom1-4/IR-limits.git
```

If you want to run numerical checks yourself, you will have to install the following dependencies:

1. [phase_space](https://github.com/tomtom1-4/phase_space) for a phase-space factorization that allows full control of the various soft and collinear parametrizations.
2. [Stripper](https://git.rwth-aachen.de/Stripper/Stripper) to extract amplitudes for various processes.
3. [Otter](Not publicly available) for Recola at quadrupole precision. Can be ignored if double-precision is sufficient. For an open source alternative, use [CutTools](https://www.ugr.es/~pittau/CutTools/) and [OneLOop](https://helac-phegas.web.cern.ch/OneLOop.html).

After installing these packages, adjust the corresponding paths in `numerics/env.sh` and source the file. Afterwards, you should be able to build the files

```bash
make
```

---

## Basic Workflow

`numerics/main` is the program that checks if the IR-limits predicted by soft currents and collinear splitting functions are actually working. It compares exact squared amplitudes obtained from `Stripper` with the approximations provided by the leading power soft or collinear factorization formula. In `numerics/main.cpp`, one has to provide a process string and the process string of the full process, i.e. including the additional radiation and some additional information like the number of resolved particles, the powers of the non-QCD couplings, the perturbative order, and a potential suffix which indicates additional model dependencies like "_QED" or "_EW". The program will then look for the process XML file inside `numerics/results/ColorCorelators`. If this file does not exist yet, one first has to run `numerics/compute_CC` with the same process. This program will then generate a Born phase-space point and will compute the relevant color and spin correlators of the process and saves the result in the XML file. This is done, so that these computations are only done once. `numerics/main` will then read this XML file to extract the phase-space point as well as the correlators. Using the `phase_space` package, the program will then construct phase-space points for the full process. **The construction of the phase-space is not fully automized**, this means that depending on the limit under consideration one has to manually select the correct clusterTree that corresponds to this limit. For nUresolved=1,2,3 this has already been implemented (see `PSF::Tree<PSF::Cluster> tree`). Depending on the limit (soft or collinear) the full phase-space point will be increasingly close to the IR limit. The program then compares the results of the IR approximations to exact results. Note that, since the computation of full phase-space points is quite computationally expensive, some results have already been saved by hand so that they are not recomputed every time. This of course only works for specific phase-space points, and can be disabled by assuring that `custom` is not set to `true`. The following IR limits were compared:
- LO soft
- NLO soft
- NNLO soft **(dipole term only)**
- LO double-soft
- NLO double-soft
- LO Triple-soft
- LO collinear
- NLO collinear
- LO triple-collinear
- NLO triple-collinear
- LO quadrupole-collinear
- (NNLO collinear not yet renormalized, and so numerical test not yet succesful)

The scale together with the relative error of the approximation are then saved into a `txt` file in `numerics/results`. The python script `numerics/results/plot.py` can be used to plot the resulting data.