# Verification record

The Project1 release branch was checked locally on 20 August 2026 using GNU
Fortran 15.2.0 on macOS.

The complete source tree compiled successfully with runtime checks enabled:

```bash
make -B \
  F90=gfortran \
  OPT="-O0 -Wall -Wextra -fcheck=all -fbacktrace" \
  EXEC=/private/tmp/sem2dsolve
```

The source tree was also compiled in release mode with `OPT="-O2"`.
`EXAMPLES/25D_TP_Dilatancy/Par.inp` completed with exit status 0 for both
builds. The test uses nonzero dilatancy, finite seismogenic width, and coupled
thermal and hydraulic diffusion. GNU Fortran reports an IEEE underflow flag at
exit for very small exponential terms; this does not stop the calculation.
