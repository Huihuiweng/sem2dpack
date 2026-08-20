# Project1 modification record

## Upstream provenance

- Original software: SEM2DPACK v2.3.8 by Jean-Paul Ampuero.
- Upstream repository: https://github.com/jpampuero/sem2dpack
- Fork lineage: https://github.com/Huihuiweng/sem2dpack
- Project1 release base: `Egao0206/sem2dpack`, branch `thermpres`, commit
  `4fca705`.

The complete Git history is retained so that the extension can be compared
with its upstream source.

## Scientific extensions

- The existing 2.5D approximation represents the finite seismogenic width in
  the elastodynamic equations.
- `SRC/bc_dynflt_tp.f90` implements one-dimensional fault-normal thermal and
  hydraulic diffusion, shear heating, thermal pressurization, and
  slip-dependent dilatancy.
- `SRC/bc_dynflt.f90` couples the thermo-hydraulic state to effective normal
  stress and fault strength, and supports temperature and pressure output.
- `SRC/bc_gen.f90` and `SRC/Makefile.depend` integrate the module with the
  SEM2DPACK boundary-condition and build systems.
- Post-processing routines include support for the additional fault outputs.

## Contribution history

The thermo-hydraulic module was adapted in part from the thermal-pressurization
implementation in MDSBI v4.1.9 by Eric M. Dunham. The module header attributes
the initial SEM2DPACK implementation to Huihui Weng and Chao Liang. The Git
history records the initial module implementation and integration by Huihui
Weng in March 2024, followed by pressure-temperature output and output-stride
changes by Yi Gao in 2024. The release retains those commits and their original
authorship metadata.

## Release preparation

The Project1 release tagged `project1-v1.0.0` adds:

- GNU Fortran portability fixes that do not alter the physical formulation;
- a compact 2.5D TP-dilatancy smoke test;
- reproducibility, provenance, citation, and build documentation;
- removal and exclusion of compiler products and runtime output files.

The production inputs used for manuscript figures belong in the associated
Open Research archive, where they can be linked to the exact software release.
