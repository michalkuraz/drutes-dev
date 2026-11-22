# NetCDF transport configuration

The active ADEnc configurations put transport and hydraulics in
`netcdf.conf`, before the concentration boundary records. Stabilization
remains in `stabilization.conf`.

## Transport and hydraulics

In `netcdf/netcdf.conf`, the section immediately after the start date sets
longitudinal dispersivity `alpha_L` [m], transverse dispersivity `alpha_T`
[m], and the first-order deposition rate `k_dep` [s⁻¹]. Use `uniform` followed by two lines: dispersivity, then deposition.
Historical unlabelled uniform records are still readable. To partition the catchment, replace
that record with a block like the following (illustrative values only):

```text
materials
200 0.2        # material 1: longitudinal/transverse dispersivity [m]
0              # material 1: deposition [1/s]
500 2.0        # material 2: dispersivity [m]
1e-5           # material 2: deposition [1/s]
end materials
```

Leave the existing Qmin record immediately after this block. As in ADE and
Richards, material IDs start at 1; provide two lines for every index through
the largest mesh material ID, including unused intermediate IDs. In a Gmsh
2.2 mesh these are the triangles’ physical surface tags, not element IDs,
boundary tags or channel path numbers. All values must be finite and
nonnegative. Omitting `alpha_T` gives isotropic dispersion; setting
`k_dep` to zero disables deposition in that material. The local suspended-mass
sink is `-k_dep * H * C` and uses the same virtual depth as storage.
Channel selection (Qmin), hydrology, and the initial concentration remain
catchment-wide.

After the channel count, add the `transport` block. Its two flags select
no-flow banks (required for ADEnc) and the accepted-step mass audit.
The conservative formulation is always used in this block. Deposition
decreases suspended mass; it is included as a sink in the mass balance.
A complete example for the current Rhine/Moselle mesh is:

```text
transport
y              # no solute flux through banks
y              # write out/adenc_mass_balance.csv
outlet profile
293200 5740000
293200 5751500
hydraulic settings
0              # lateral sources: 0 = none, 1 = explicit lateral.conf
1              # forcing: 0 = daily steps, 1 = linear interpolation
1e-10          # relative water-balance tolerance
20000          # maximum projection PCG iterations
1.0            # maximum relative edge-flux correction
y              # remove whole disconnected components with neither inlet nor outlet
end transport
```

Each point is `x y` in the mesh CRS (metres). Order them so that downstream
is to the **left** when looking from the first point to the second. ADEnc
looks for midpoints of open active-domain edges in a band extending 2.5%
of the segment length on either side. The routed flow must point downstream
and outward across the edge. Place the segment along the outlet face,
slightly overlapping the real opening. This selects existing edges,
rather than cutting or retriangulating the mesh. Where two eligible edges
project onto the same part of the profile, the closer edge takes priority.
Reversed or misplaced points can fail initialization.

The tagged edges on the exterior of the mesh supply inlets. Their discharge
comes from `Qrouted`, with the existing length weighting per inlet ID.
Outlets use the existing balanced-discharge projection. No node IDs are
entered. Missing elevation is not used to choose ports.

The `hydraulic settings` subsection is optional. Defaults are
`0, 1, 1e-10, 20000, 1.0, n` (in that order). ADEnc's conservative
formulation requires no-flow banks and standard 2D Picard.
Omit `outlet profile` and its settings entirely to use unreconstructed
flux with a warning. A supplied profile which cannot be resolved is an
error, not a silent fallback. The independent mRM routing input is still
required in both modes.

The root configuration and scenarios 1, 2, and 3 keep their previous
unreconstructed transport. Scenarios 3a and 3b and the runnable examples
specify the profile above.

## Stabilization

`stabilization.conf` has exactly four data records:

```text
# Enable SUPG
y
# SUPG multiplier
2.0
# Enable residual viscosity (shock capturing)
y
# Residual-viscosity multiplier
1.0
```

| Variant | SUPG flag | SUPG multiplier | Residual-viscosity flag | Multiplier |
|---|---|---|---|---|
| Galerkin | n | 2.0 | n | 1.0 |
| SUPG | y | 1.0 | n | 1.0 |
| SUPG2 | y | 2.0 | n | 1.0 |
| SUPG2 + residual viscosity | y | 2.0 | y | 1.0 |

Both multipliers must be finite and nonnegative, even when disabled.
Zero disables the corresponding contribution. The two switches are independent.
An absent file disables both stabilizations with multiplier defaults of one.
Neither option guarantees positivity; residual viscosity can increase spreading.
The physical dispersivity values remain in `netcdf.conf`.

Flags must be lowercase `y` or `n` on separate records. Blank lines and `#`
comments are allowed. Changing these files does not require recompilation.

## Migrating an existing configuration

Keep the no-flow banks and mass-audit flags in the block after the channel count.
Remove the old conservative-transport flag; set the bank flag to `y`.
If reconstruction was used, replace the node-pair list with an outlet
profile and keep any non-default hydraulic settings. Check the physical
inlet tags in the mesh: they must mark only the actual inlet openings.
Archive the old `transport.conf` outside the configuration; a merged block
and the old file cannot coexist. Unmerged historical inputs remain readable
for reproducibility.

Every new model run needs its own empty `out/`: main clears `out/*` at startup.
Numerical calculations, boundary formulas, flux projection and stabilization
formulas were not changed by this input consolidation.
