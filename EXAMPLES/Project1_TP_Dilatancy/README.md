# Project1 TP-dilatancy smoke test

This compact example checks the coupled 2.5D dynamic-fault implementation with
thermal pressurization, fault-normal diffusion, and nonzero dilatancy. It is
intended to verify a build, not to reproduce a manuscript figure.

Build `SRC/sem2dsolve`, then run it from this directory. The test normally
finishes in a few seconds to tens of seconds on a desktop system and should
exit with status 0.

The checking-phase PostScript diagnostic is disabled through the third entry
of `verbose` because that legacy plotting path is not portable to all current
systems.
