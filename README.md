# SEM2DPACK TP-dilatancy extension

This repository is a research fork of SEM2DPACK v2.3.8. It extends the
dynamic-fault implementation with a finite-width 2.5D approximation and a
fault-normal thermo-hydraulic module that includes thermal pressurization,
hydraulic diffusion, and shear-induced dilatancy.

The `thermpres` branch contains the implementation used for the dynamic rupture
simulations associated with the Project1 manuscript. The corresponding
software release is identified by the tag `project1-v1.0.0`. Production input
files and processed outputs from the full parameter sweeps are not included.

## Main extension

The `BC_DYNFLT_TP` input block activates the module in
`SRC/bc_dynflt_tp.f90`. The module evolves temperature and pore-pressure
perturbations on a one-dimensional grid normal to the fault. Shear heating
provides the thermal source, while a slip-dependent inelastic porosity change
provides the dilatancy source. The ratio `Phi/beta` controls the pressure scale
of dilatant suction, and all input quantities must use a consistent unit
system.

Set `tp_file=.true.` to write the full fault-normal temperature and pressure
fields. Two representative input examples are included:

- TP-only: `EXAMPLES/Thermpres_SWF/Par.inp`
- TP-dilatancy: `EXAMPLES/25D_TP_Dilatancy/Par.inp`

The TP-dilatancy example provides a small 2.5D smoke test with both processes
enabled.

## Build

The release has been compiled with GNU Fortran 15.2.0. From the repository
root:

```bash
cd SRC
make F90=gfortran OPT="-O3" EXEC="$PWD/sem2dsolve"
```

The original Makefile also supports Intel Fortran through its user settings.

## Run the smoke test

```bash
cd EXAMPLES/25D_TP_Dilatancy
../../SRC/sem2dsolve
```

Successful completion ends with `Program SEM2DPACK: end`. This compact test
checks the 2.5D elastic model, slip-weakening and time-weakening fault laws,
and the coupled TP-dilatancy update. It is not a production manuscript model.

## Provenance and license

The original SEM2DPACK software was developed by Jean-Paul Ampuero and is
distributed under the GNU General Public License, version 2 or later. The
original `Copyright`, `LicenseNotice`, and `GeneralPublicLicense` files are
retained unchanged. Parts of the thermo-hydraulic implementation were adapted
from MDSBI v4.1.9 by Eric M. Dunham. See `MODIFICATIONS.md` for the extension
history and file-level attribution.

Original software citation:

J. P. Ampuero (2012), SEM2DPACK, a spectral element software for 2D seismic
wave propagation and earthquake source dynamics, v2.3.8. Zenodo.
https://doi.org/10.5281/zenodo.230363
