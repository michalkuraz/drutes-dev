"""Bounds-checked tests of the production containment function; never run main."""
from pathlib import Path
import shlex
import shutil
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def run(command, cwd):
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=300)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


@pytest.fixture(scope='module')
def driver(tmp_path_factory):
    needed = {name: shutil.which(name) for name in ('gfortran', 'make', 'nf-config', 'nc-config')}
    if not all(needed.values()):
        pytest.skip('Fortran and NetCDF development tools required for production-module tests')
    work = tmp_path_factory.mktemp('inside-production')
    build = work / 'build'
    nf_flags = shlex.split(run([needed['nf-config'], '--fflags'], work))
    flags = shlex.join(['-DHAVE_NETCDF', '-fimplicit-none', '-fcoarray=single', '-fdefault-real-8',
                        '-fcheck=all', '-fbacktrace', '-O0', '-g', f'-J{build}/mods', *nf_flags])
    run([needed['make'], 'build_target', f'BUILD={build}', f'BINDIR={work}/bin', f'FFLAGS={flags}'], ROOT)
    objects = sorted(p for p in (build / 'objs').glob('*.o') if p.name != 'main.o')
    libraries = [f"-L{run([needed['nc-config'], '--libdir'], work).strip()}",
                 *shlex.split(run([needed['nf-config'], '--flibs'], work))]
    executable = work / 'test-inside'
    run([needed['gfortran'], '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
         '-ffpe-trap=invalid,zero,overflow', '-I', str(build / 'mods'), *nf_flags,
         str(ROOT / 'tests/fortran/test_geom_inside.f90'), *map(str, objects),
         *libraries, '-o', str(executable)], work)
    return executable, work


def test_random_geometry_and_boundaries(driver):
    executable, work = driver
    assert 'Deterministic inside checks passed' in run([str(executable)], work)


def test_unique_random_rhine_triangles(driver):
    executable, work = driver
    mesh = ROOT / 'drutes.conf1/mesh/mesh.msh'
    if not mesh.exists():
        pytest.skip('Local Rhine mesh fixture not present')
    # Preserve input triangle order, matching read_2dmesh_gmsh's FE indices.
    rows = mesh.read_text().splitlines()
    start = rows.index('$Nodes')
    nodes = {int(v[0]): [float(v[1]), float(v[2])]
             for v in (s.split() for s in rows[start+2:start+2+int(rows[start+1])])}
    start = rows.index('$Elements')
    ids = [list(map(int, v[-3:])) for v in
           (s.split() for s in rows[start+2:start+2+int(rows[start+1])]) if v[1] == '2']
    triangles = np.array([[nodes[n] for n in element] for element in ids])
    rng = np.random.default_rng(1729)
    selected = rng.integers(len(triangles), size=1000)
    weights = rng.dirichlet([3, 3, 3], size=len(selected))
    # Generate strict interior points from known elements, not from inside itself.
    vertices = triangles[selected]
    queries = vertices[:, 0] + weights[:, 1, None]*(vertices[:, 1]-vertices[:, 0]) + \
              weights[:, 2, None]*(vertices[:, 2]-vertices[:, 0])
    known = [(10679, 405591.710404, 5352097.992453),
             (6415, 461843.491141, 5510644.713862),
             (10417, 361575.927733, 5627442.674994),
             (17993, 329703.451993, 5725439.007906)]
    fixture = work / 'rhine-queries.txt'
    with fixture.open('w') as stream:
        stream.write(f'{len(triangles)} {len(selected)+len(known)}\n')
        np.savetxt(stream, triangles.reshape(-1, 6), fmt='%.17e')
        np.savetxt(stream, np.column_stack((selected+1, queries)), fmt=['%d', '%.17e', '%.17e'])
        for element, x, y in known:
            stream.write(f'{element} {x:.17e} {y:.17e}\n')
    output = run([str(executable), str(fixture)], work)
    assert 'Rhine unique triangle queries passed:' in output
