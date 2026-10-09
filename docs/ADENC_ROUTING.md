# Mandatory mHM routing network in ADEnc

Every ADEnc initialization now requires the fixed input:

```
drutes.conf/netcdf/mRM_restart_001.nc
```

Missing, unreadable, incompatible or invalid networks cause `ERROR STOP`.
There is no enable switch or configurable routing filename. This requirement
is specific to ADEnc; other DRUtES models are unchanged. Existing ADEnc input
sets must supply this file before using the new executable. A symlink is valid
and recommended when several cases share one mHM domain. Existing saved
results must not be overwritten/recomputed merely to add this input.

## Filename convention and provenance

mHM does **not** universally require the suffix `_001`. Its namelist explicitly
sets `mrm_file_RestartOut(domain)` and `mrm_file_RestartIn(domain)`. The documented
example names domain 1 `mRM_restart_001.nc`, domain 2 `mRM_restart_002.nc`.
ADEnc currently supports one routing domain and intentionally adopts domain
1's conventional filename. Supply the corresponding domain's data, not an
unrelated file renamed to satisfy the existence check.

References:

- https://mhm.pages.ufz.de/mhm/latest/mhm_8nml_source.html
- https://mhm.pages.ufz.de/mhm/latest/mo__common__read__config_8_f90_source.html
- https://mhm.pages.ufz.de/mhm/stable/namespacemo__mrm__net__startup.html

The Rhine source on hydrocalc is:
`/mnt/stock/maryam/rhein/restart/mRM_restart_001.nc`.
Its morphological files are `../morph/{dem,fdir,facc}.asc` and its routing
resolution is 0.125 degrees. The network has 1439 nodes and 1438 links.

During implementation a 125 KiB **routing-only extract** was downloaded before
VPN/DNS access failed. The local fixed-name input contains these exact copied
static routing fields, not the complete 467 MiB mHM restart or its evolving
water states. Its `source_file` attribute identifies the original. It must
not be used to restart mHM itself. ADEnc only needs static routing variables;
the reader accepts either this extract or the full original restart. Once
server access returns, the full original can replace the local extract without
changing the interface. NetCDF data files are excluded from Git; distribute
the input separately and install/symlink it into each case.

## Implementation/API (`src/models/fluxLS/ncrouting.f90`)

- `read_adenc_routing()` loads the mandatory fixed path, validates it and logs
  node/outlet/confluence counts.
- `routing_read(path, ok, message)` is the lower-level reader for tests/tools.
  Failure clears the cache. It does not restore mHM dynamic restart states.
- `routing_cell(cell, downstream_cell, upstream_cells, direction, fdir, ok)`
  returns all upstream cells, the unique downstream cell and a unit UTM
  direction between cell centres. `cell` and neighbour IDs are DRUtES
  hydrological **array indices**, not mHM IDs, Gmsh tags or FE indices.
- `routing_element(element, ...)` queries the same data via
  `el2ncgrid(element)`. `element` is the FE array index.
- `direction = routing_flow_direction(element, ok)` is a convenience function
  returning just the two-component unit UTM vector. Always inspect `ok`.
- `routing_apply_directions(ok, message)` installs those vectors in
  `ncfluxdata%fluxvct` after channel reading and FE mapping, **before** width
  preparation and hydraulic reconstruction. `routing_close()` clears the cache.

The reader uses `L11_Id`, `L11_domain_mask`, `L11_domain_lat`, `L11_domain_lon`,
`L11_fDir`, `L11_fromN`, `L11_toN`. It validates variable ranks and common grid
dimension IDs, masks, dense unique mHM IDs, coordinate matching to Qrouted
centres (absolute tolerance 1e-7 degrees), valid link endpoints, neighbouring
cells, one outgoing link per non-outlet, and absence of cycles. A trailing
paired fill-value link, as present in the Rhine restart, is ignored. Outlets
have `L11_fDir=0` and no outgoing link.

Mapping uses **actual latitude/longitude values**, not the restart's misleading
dimension names, flattening assumptions, or matching mHM node IDs with FE IDs.
Descending/reordered axes are supported. Query outputs use
`cell=(ilat-1)*nlon+ilon`, matching the hydrological quadrilateral numbering.
The raw `fdir` code is exposed for diagnostics but not interpreted as a generic
GIS D8 encoding. Explicit `fromN -> toN` links determine physical directions.

At a true routing outlet the query succeeds with downstream 0 and direction
0. For actual FE flow at that outlet, the existing explicitly ordered channel
polyline supplies the exit direction; this exception is logged. Unmapped or
masked fringe FE cells retain existing zero-width handling. No missing
downstream link is silently treated as an outlet.

## Scope and scientific limits

Network directions replace nearest-polyline directions for mapped non-outlet
FE cells. This can change effective widths, storage and the raw field used by
RT0 reconstruction; old concentration outputs were computed with the old
directions. Adding a routing file to an archived case does **not** retroactively
validate those results or make them outputs of this implementation.

This does not infer lateral rates from Q differences, restore historical
storage, transfer upstream solute across unresolved tributaries, automatically
open FE boundary ports, or fix the case-2 rising-wave storage deficit. The
existing explicit lateral-source and port formulations remain separate.
No generic FEM, main solver or Schwarz numerical code was changed by routing.

## Tests

`python -m pytest -q tests/test_ncrouting.py` compiles the production reader
with bounds checking and tests reordered grids/IDs, branching, FE queries,
directions, cache resets, mandatory missing inputs, malformed variables,
invalid coordinates/IDs/masks/codes/links and cycles. Tests never call main.

`tests/fortran/test_ncrouting.f90 <routing-file> real <Qrouted-file>` checks
all real links/upstream reciprocity and writes `routing-analysis.csv` in its
working directory; run it in a **new diagnostic directory**. It never starts
a transport simulation. Full initialization/assembly checks must likewise
use new directories because output routines create files.
