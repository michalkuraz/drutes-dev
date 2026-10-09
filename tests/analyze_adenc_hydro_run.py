"""Read-only screening and independent P0(H)*P1(C) inventory of short hydro tests.

Uses standard Python only. Diagnostics must be exported by
export_adenc_hydro_domain.f90 from the same mesh/date/filtered geometry.
No main invocation, model modification, clipping or conservation claim beyond
the supplied numerical field/boundary budget. Within the first forcing day only.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def nodal_blocks(path):
    """Read all scalar Gmsh2 NodeData blocks, including multiple times per file."""
    with path.open() as stream:
        for line in stream:
            if line.strip() != '$NodeData':
                continue
            strings = [next(stream).strip() for _ in range(int(next(stream)))]
            reals = [float(next(stream)) for _ in range(int(next(stream)))]
            ints = [int(next(stream)) for _ in range(int(next(stream)))]
            if len(ints) < 3 or ints[1] != 1 or not reals:
                raise ValueError('Expected scalar NodeData with time and record count')
            values = {}
            for _ in range(ints[2]):
                row = next(stream).split()
                if len(row) != 2 or int(row[0]) in values:
                    raise ValueError('Invalid or duplicate scalar nodal record')
                values[int(row[0])] = float(row[1])
            if next(stream).strip() != '$EndNodeData':
                raise ValueError('Missing EndNodeData')
            yield reals[0], strings, values


def analyze(folder, geometry, constant=None):
    folder, geometry = Path(folder), Path(geometry)
    nodes = {}
    for row in (geometry/'nodes.txt').read_text().splitlines():
        v = row.split(); nodes[int(v[0])] = (float(v[1]), float(v[2]))
    next_q = {}
    for row in (geometry/'next.txt').read_text().splitlines():
        v = row.split(); next_q[int(v[0])] = float(v[1])
    triangles, active_nodes = [], set()
    for row in (geometry/'active.txt').read_text().splitlines():
        v = row.split(); e = int(v[0]); ids = tuple(map(int, v[1:4]))
        xyz = [nodes[i] for i in ids]
        area = abs((xyz[1][0]-xyz[0][0])*(xyz[2][1]-xyz[0][1])-
                   (xyz[1][1]-xyz[0][1])*(xyz[2][0]-xyz[0][0]))/2
        width, q = float(v[6]), float(v[7]); q1 = next_q[e]
        # The staged Rhine constants are vref1.7m/s, Qref1672m3/s.
        h0 = q/(width*1.7*(q/1672)**.4)
        h1 = q1/(width*1.7*(q1/1672)**.4)
        triangles.append((ids, area, h0, h1)); active_nodes.update(ids)
    with (folder/'out/adenc_mass_balance.csv').open() as stream:
        budget = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]
    snapshots = []
    files = sorted((folder/'out').glob('ADE_in_watershed_solute_concentration-*.msh'))
    files = [f for f in files if f.stem.rsplit('-', 1)[1].isdigit() and '-el_avg-' not in f.name]
    for path in files:
        for time, names, values in nodal_blocks(path):
            if not 0 <= time <= 86400:
                raise ValueError('This diagnostic supports only the first forcing day')
            if not active_nodes <= values.keys():
                raise ValueError('Missing active FE nodes')
            c = [values[i] for i in active_nodes]
            mass = math.fsum(area*(h0+(h1-h0)*time/86400)*
                             math.fsum(values[i] for i in ids)/3
                             for ids, area, h0, h1 in triangles)
            negatives = math.fsum(area*(h0+(h1-h0)*time/86400)*
                                  math.fsum(max(-values[i], 0) for i in ids)/3
                                  for ids, area, h0, h1 in triangles)
            entry = dict(time_s=time, active_min=min(c), active_max=max(c),
                         finite=all(math.isfinite(v) for v in c), inventory=mass,
                         negative_nodal_inventory=negatives,
                         nodes_above_1p001=sum(v>1.001 for v in c),
                         nodes_below_minus1e_3=sum(v<-.001 for v in c),
                         positive_nodes=sum(v>1e-6 for v in c))
            if constant is not None:
                entry['max_constant_error'] = max(abs(v-constant) for v in c)
            matches = [r for r in budget if abs(r['time_s']-time)<1e-7]
            if matches:
                entry['spatial_minus_audit_inventory'] = mass-matches[-1]['inventory']
            snapshots.append(entry)
    snapshots.sort(key=lambda r: r['time_s'])
    terminal = (folder/'terminal.log').read_text(errors='replace')
    return dict(folder=str(folder), exit_status=int((folder/'exit.status').read_text()),
                numerical_finish='F I N I S H E D' in terminal,
                accepted_steps=len(budget), last_time_s=budget[-1]['time_s'],
                max_relative_budget_error=max(abs(r['relative_error']) for r in budget),
                max_free_residual=max(abs(r['free_residual']) for r in budget),
                initial_inventory=budget[0]['initial_inventory'],
                final_budget=budget[-1], active_elements=len(triangles),
                active_nodes=len(active_nodes), snapshots=snapshots,
                notes=['Short numerical check, not calibrated Rhine physics or positivity validation.',
                       'No lateral sources; deliberately restricted test outlet.',
                       'P0 centroid-filtered fixed geometry; linear daily H interpolation.',
                       'Independent spatial inventory uses saved P1 nodal C and reconstructed P0 H.',
                       'Dirichlet inventory budget includes dispersive and algebraic boundary transfer; not just Q*C.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder'); parser.add_argument('geometry')
    parser.add_argument('--constant', type=float)
    args = parser.parse_args()
    print(json.dumps(analyze(args.folder, args.geometry, args.constant), indent=2, allow_nan=False))
