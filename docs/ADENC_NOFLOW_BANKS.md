# ADEnc impermeable solute banks

Implemented 2026-10-08. Enable with `drutes.conf/netcdf/riverbank.conf` containing
`y`. Root configuration is enabled. A missing file or `n` retains the legacy
absorbing C=0 mask. Existing run/export configuration copies are unchanged.
Supported: ADEnc, 2D P1 triangles, standard Picard (`it_method=0`). Schwarz
modules still compile/link, but this feature rejects Schwarz execution.

## Domain and boundaries

The initial Q/Qmin node classification still defines the same fixed active
element mask (at least one qualifying node per triangle). This is a FEM-grid
approximation, not a conforming remesh at exact hydrological cell borders.
Only active triangles assemble capacity, transport, and source terms.
All nodes participating in these triangles regain unknown concentration DOFs,
except nodes with original physical Dirichlet boundary IDs. Nodes belonging
only to unused elements remain eliminated using the reserved `addedbc` ID;
that ID is forced to constant zero solely for their output/DOF handling.

Undirected edge hashing identifies active/inactive interfaces. Shared
active/active edges are never banks. The true exterior of the original FEM
mesh is distinct and keeps its existing physical/natural boundary treatment;
it is not automatically sealed by this feature. Matching original boundary
IDs identify port edges and exclude them from bank treatment. Original
inlet/outlet boundary callbacks and data are unchanged.

Legacy `terrain_slopes` also assigned addedbc to triangles with missing DEM.
The participating-node restoration is performed after that annotation: missing
DEM no longer silently imposes C=0 in this mode. ADEnc convection directions
are subsequently assigned from channel polylines. Existing DEM coverage
warnings remain; terrain geometry has not been repaired by this change.

Nodal hydrological output now selects an adjacent active triangle if one
exists, rather than blindly selecting the first (possibly inactive) triangle.
This does not resolve all previously noted observation-history discrepancies.

## Zero total solute flux

On an internal bank with outward active-domain normal n, prescribe

    J.n = (q*C - K*grad(C)).n = 0.

With the current advective volume form H*C_t + q.grad(C) = div(K grad(C)),
this becomes a Robin condition K*grad(C).n = (q.n)*C. DRUtES uses negative
capacity/stiffness signs. Its edge addition is therefore

    +dt * integral_edge (q.n) Ni Nj ds.

Two-point Gauss integration evaluates the edge matrix. Hydrological Q is
queried on the active-element side (a 1e-8 barycentric inward offset avoids
ambiguous exact-cell-boundary lookup); q=Q/W times the stored channel direction.
Width and daily time lookup match the volume callback conventions. Failed/dry
lookups and zero-width elements have zero advective edge contribution, matching
the existing zero-convection behavior. For tangent flow this addition vanishes:
zero normal dispersion is a natural condition. No Newton solver is needed.

The bank correction is applied after the existing SUPG/shock element hook.
It is the natural total-flux condition for the effective diffusion operator
including residual artificial diffusion, when that is enabled.

## Verification

- Optimized project build succeeds, including both Schwarz modules.
- Full pytest suite: 35 passed (2026-10-08).
- Production-module bounds-checked Fortran test checks edge integrals/signs,
  tangent flow, active/active exclusion, active/inactive vs external edge
  classification, port preservation, unknown bank DOFs, and inactive-first
  nodal adjacency.
- A manufactured closed square uses the production assembler and bank hook,
  synthetic cached constant hydrology and anisotropic diffusion. A separate
  inactive triangle raises an error if accidentally assembled. Both lumped
  and consistent Euler capacity, Galerkin/SUPG2/SUPG2+shock1 are exercised for
  30 steps each (180 assembled steps). Constant-H inventory stays within 1e-11;
  the global column-sum identity is checked to 1e-12. The test-only option
  `close_exterior=.true.` closes the manufactured external box; production
  does not use that option.
- Rhine production initialization, staged from the exported SUPG2 case with
  only riverbank.conf added, succeeds without invoking main or time stepping:
  4,968 active elements, 2,869 participating nodes, 723 internal bank edges,
  15,099 unused nodes, 46 inlet nodes, inlet series shape (3,2).
  There are 713 active FE elements with zero effective width, inherited from
  the existing approximate mask/width mapping. They are not newly repaired.

## Remaining limitations / before a new Rhine run

This change removes the absorbing C=0 internal banks, not all conservation
issues. H_t+div(q), discontinuous interior flow traces, daily storage changes
and inlet/outlet budgets still require a discrete mass audit. A Robin condition
can cancel a nonzero normal water flux in the solute equation; it does not make
the hydrological vector field tangential to the raster/FEM bank.

In the tested first-scenario mesh only physical inlet ID101 is present;
ID102 is the RESERVED UNUSED ID, not an explicit outlet. An active-domain
termination inside the full mesh will therefore be a sealed internal bank.
If that termination should be an open river outlet, it must be explicitly
identified/labelled and supplied with an appropriate outlet condition before
interpreting a long run physically. Do not relabel every bank as an outlet.
Natural zero dispersive flux on an unlabelled external mesh edge is not the
same as this internal-bank condition.

No production time simulation or server deployment was performed by this
implementation. Current server executables/configurations/results are intact.
