"""Read-only area audit of complete scalar Gmsh 2.2 NodeData snapshots.

Integrates the P1 concentration, NOT contaminant mass. For variable storage H,
mass is integral(H*C); a boundary-flux budget is also required for conservation.
Negative/positive integrals are exact for the represented affine triangle field
(clip each triangle at C=0); no concentration clipping changes the solution.
Usage: python3 integrate_adenc_concentration.py CASE_DIR [CASE_DIR ...]
Incomplete snapshots in an actively written file are ignored.
"""
import argparse
import csv
import math
from pathlib import Path
import sys


def mesh(rows):
    i = rows.index('$Nodes')
    count = int(rows[i + 1])
    nodes = {}
    for row in rows[i + 2:i + 2 + count]:
        r = row.split()
        nodes[int(r[0])] = (float(r[1]), float(r[2]))
    i = rows.index('$Elements')
    count = int(rows[i + 1])
    triangles = []
    for row in rows[i + 2:i + 2 + count]:
        r = list(map(int, row.split()))
        if r[1] == 2:
            ids = r[3 + r[2]:]
            if len(ids) != 3:
                raise ValueError('Only three-node triangles are supported')
            triangles.append(tuple(ids))
    if not triangles:
        raise ValueError('No P1 triangles found')
    return nodes, triangles


def snapshots(rows):
    i = 0
    while i < len(rows):
        if rows[i] != '$NodeData':
            i += 1
            continue
        try:
            j = i + 1
            ns = int(rows[j]); j += 1 + ns
            nr = int(rows[j]); j += 1
            real_tags = list(map(float, rows[j:j + nr])); j += nr
            ni = int(rows[j]); j += 1
            tags = list(map(int, rows[j:j + ni])); j += ni
            if len(tags) < 3 or tags[1] != 1:
                raise ValueError('Expected scalar NodeData')
            n = tags[2]
            if j + n >= len(rows) or rows[j + n] != '$EndNodeData':
                break
            values = {}
            for row in rows[j:j + n]:
                r = row.split()
                values[int(r[0])] = float(r[1])
            if len(values) != n or not all(math.isfinite(v) for v in values.values()):
                raise ValueError('Duplicate node IDs or nonfinite values')
            yield real_tags[0], values
            i = j + n + 1
        except IndexError:
            break


def area(a, b, c):
    return abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])) / 2


def positive_integral(vertices):
    """Exact integral of max(C,0) on a triangle, by zero-line clipping."""
    polygon = []
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        if a[2] >= 0:
            polygon.append(a)
        if (a[2] < 0 < b[2]) or (b[2] < 0 < a[2]):
            f = a[2] / (a[2] - b[2])
            polygon.append((a[0]+f*(b[0]-a[0]), a[1]+f*(b[1]-a[1]), 0.0))
    return math.fsum(area(polygon[0], polygon[k], polygon[k+1]) *
                     (polygon[0][2]+polygon[k][2]+polygon[k+1][2])/3
                     for k in range(1, len(polygon)-1))


def audit(case):
    path = case / 'out/ADE_in_watershed_solute_concentration-0.msh'
    rows = path.read_text().splitlines()
    nodes, triangles = mesh(rows)
    areas = [area(*(nodes[n] for n in tri)) for tri in triangles]
    if any(a <= 0 for a in areas):
        raise ValueError('Degenerate triangle')
    total_area = math.fsum(areas)
    for time, values in snapshots(rows):
        signed, pos, neg = [], [], []
        for tri, a in zip(triangles, areas):
            v = [(*nodes[n], values[n]) for n in tri]
            signed.append(a * math.fsum(p[2] for p in v) / 3)
            pos.append(positive_integral(v))
            neg.append(positive_integral([(x, y, -c) for x, y, c in v]))
        s, p, n = map(math.fsum, (signed, pos, neg))
        if not math.isclose(s, p-n, rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError('Signed/positive/negative integrals disagree')
        yield [case.name, time, time/86400, total_area, s, p, n,
               n/p if p else 0, min(values.values()), max(values.values())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases', nargs='+', type=Path)
    args = parser.parse_args()
    writer = csv.writer(sys.stdout)
    writer.writerow(['case', 'time_s', 'time_days', 'area_m2', 'integral_C_dA',
                     'integral_positive_C_dA', 'integral_negative_magnitude_dA',
                     'negative_positive_ratio', 'Cmin', 'Cmax'])
    for case in args.cases:
        writer.writerows(audit(case))


if __name__ == '__main__':
    main()
