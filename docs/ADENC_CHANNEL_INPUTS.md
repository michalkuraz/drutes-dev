# ADEnc channel paths and inlet boundaries

Updated 2026-10-01.

In `drutes.conf/netcdf/netcdf.conf`, the integer after the initial concentration
is the number of channel paths. A value of 1 loads `channel.dat`; larger
values additionally load `channel2.dat`, `channel3.dat`, through `channelN.dat`
from the same directory.

Each file contains at least two points with exactly two whitespace-separated
finite numbers per data line: projected x and y, in the mesh coordinate system.
Spaces and tabs, blank lines, and comments introduced by `#` are supported.
Points are ordered downstream: consecutive points determine the flow direction.
Duplicate consecutive points within one file are rejected. The same coordinate
in different files is allowed, for example at a confluence.

Separate files contribute separate segments. No segment is created from the
end of one file to the beginning of another. The existing nearest-segment
direction rule considers all loaded segments; intersecting or closely spaced
branches still require checking the assigned directions around the confluence.

Channel count and inlet-boundary count are independent. Positive mesh boundary
IDs start at 101. ADEnc reserves `max(original mesh boundary ID)+1` for inactive
nodes, before the hydrological/DEM masks are applied. The configuration must
provide one boundary record for every ID from 101 through that reserved ID,
in ascending order.

For two physical inlets labeled 101 and 102, the records are:

```text
101 1 y 0.0
102 1 y 0.0
103 1 n 0.0
```

The first two records load independent `101.bc` and `102.bc` concentration
series from `drutes.conf/netcdf/`. The last record sets the inactive-node value.
The reserved boundary record is required even if all mesh nodes are active.
Time-series files follow the existing DRUtES reader contract; include a final
time strictly beyond all possible boundary evaluations (including the last
step) because the current Dirichlet callback selects the preceding data row.

Existing one-inlet configurations retain records 101 and 102. This code extension
does not itself add a tributary polyline, change the FE mesh, or assign a new
physical inlet to a prepared scenario.
