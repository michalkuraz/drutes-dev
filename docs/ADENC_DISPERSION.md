# Direction-aligned ADEnc dispersion

Implemented 2026-10-07. The dispersivity record in
`drutes.conf/netcdf/netcdf.conf` accepts either:

```
200.0 0.2  # alpha_L [m], alpha_T [m]
```

or one legacy scalar, which preserves isotropic dispersion (`alpha_T=alpha_L`).
Both values must be finite and nonnegative, on the same physical line.
The following Qmin, initial concentration, channel count and boundary records
keep their existing positions. Initialization logs both dispersivities.

## Formulation and dimensions

The existing convection vector is the depth-integrated water flux
`q = (Q/W) d`, with unit channel direction `d`, discharge Q [m3/s] and effective
width W [m]. The dispersion callback now returns the symmetric tensor

`K = |q| [alpha_T I + (alpha_L-alpha_T) d d^T]` [m3/s].

Off-diagonal entries are retained. Zero flux gives a zero tensor. Equal
dispersivities recover the old isotropic tensor. The scalar callback reports
the longitudinal coefficient `alpha_L |q|`; it is not a replacement for K.
ADEnc remains a 2D integration-point-based callback.

The existing storage coefficient is effective depth `h=Q/(W v)` [m]. Thus
the local diffusion coefficients after division by h are `D_L=alpha_L v`
and `D_T=alpha_T v` [m2/s]. Do not replace q by velocity inside K without also
changing the depth-integrated governing equation. Storage and convection
were not changed here; no molecular diffusion was added.

The directional tensor structure follows the mechanical-dispersion formulation
described in the [USGS MODFLOW 6 transport documentation](https://pubs.usgs.gov/tm/06/a61/tm6a61.pdf).
That reference supports the tensor algebra, not river-specific calibration.

## Proposed starting values, not calibration

Root defaults are alpha_L=200 m and alpha_T=0.2 m. At a hypothetical speed
of 1 m/s these imply 200 and 0.2 m2/s; the actual model speed is variable.
This deliberately reduces along-channel spreading tenfold relative to the
previous 2000 m value and strongly limits cross-channel spreading.
These are exploratory effective parameters, not inferred from the finished run
and not a validated parameterization of the Rhine or its coarse grid-cell width.
Longitudinal and transverse mixing can differ by orders of magnitude:
[USGS Missouri tracer measurements](https://pubs.usgs.gov/publication/wsp1899G)
report approximately 1490 and 0.12 m2/s respectively. Those measured values
must not be transferred directly to the Rhine.

Follow-up should separate sensitivity to alpha_L (e.g. 100, 200, 500 m) from
alpha_T (e.g. 0.1, 0.2, 1 m), inspect time-step/mesh convergence, concentration
bounds and mass balance, and compare with tracer data if available. Lower
dispersion can expose convection-dominated FEM oscillations. It does not fix
the observed mismatch between point 1's exported history and spatial solution.

Only root defaults and the root executable were updated. Existing scenario
inputs `drutes.conf1`, `drutes.conf2`, `drutes.conf3`, `drutes.conf3a`,
`drutes.conf3b`, finished run inputs/outputs, and GUI project copies remain
unchanged for reproducibility. To use anisotropy in a new staged scenario,
set its dispersivity record explicitly to `200.0 0.2` and copy the new binary.
No transport simulation was launched as part of this change.

## Verification

`make build_target` succeeded with NetCDF support. The complete test suite
passed (35 tests). New checks cover single/pair record parsing, invalid and
nonfinite inputs, preservation of the following Qmin record, eigenvectors and
eigenvalues for rotated/reversed flow, symmetry, the legacy isotropic case,
zero flux and zero dispersivity, and the real `ADElsdisp` PDE callback linked
against the production modules with bounds checking. FEM stiffness assembly
already uses full matrix multiplication with K, including off-diagonal terms.
These are implementation checks, not transient numerical validation.
