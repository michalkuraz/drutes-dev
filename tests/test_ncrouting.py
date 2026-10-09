"""Production mRM reader, reordered grids/IDs, confluences and invalid inputs."""
from pathlib import Path
import shlex
import shutil
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def fixture(path, change=None):
    nc = pytest.importorskip('netCDF4')
    # Deliberately different order from Qrouted; cannot equate IDs with FE/cell IDs.
    ids = np.array([[3, 4], [2, 1]], dtype='i4')
    fdir = np.array([[0, 4], [1, 4]], dtype='i4')
    lat = np.array([[50.875, 50.875], [51, 51]], dtype='f8')
    lon = np.array([[7.125, 7], [7.125, 7]], dtype='f8')
    mask = np.ones((2, 2), dtype='i4')
    links = np.array([4, 1, 2, -9999], dtype='i4')
    targets = np.array([3, 2, 3, -9999], dtype='i4')
    if change == 'cycle':
        targets[2] = 1
    elif change == 'duplicate-link':
        links[0] = 1
    elif change == 'bad-target':
        targets[0] = 99
    elif change == 'missing-link':
        links[2] = targets[2] = -9999
    elif change == 'self-loop':
        targets[0] = 4
    elif change == 'outlet-link':
        fdir[0, 0] = 4
    elif change == 'duplicate-id':
        ids[0, 1] = 1
    elif change == 'invalid-id':
        ids[0, 1] = 99
    elif change == 'invalid-mask':
        mask[0, 0] = -9999
    elif change == 'invalid-code':
        fdir[0, 1] = 3
    elif change == 'grid-mismatch':
        lat += .01
    elif change == 'nan':
        lon[0, 0] = np.nan
    with nc.Dataset(path, 'w') as ds:
        ds.createDimension('rows', 2)
        ds.createDimension('cols', 2)
        ds.createDimension('links', 4)
        for name, values in [('L11_Id', ids), ('L11_domain_mask', mask),
                             ('L11_fDir', fdir), ('L11_domain_lat', lat), ('L11_domain_lon', lon)]:
            if change == 'missing-variable' and name == 'L11_fDir':
                continue
            dims = ('cols', 'rows') if change == 'swapped-variable-dims' and name == 'L11_fDir' else ('rows', 'cols')
            var = ds.createVariable(name, values.dtype, dims, fill_value=-9999)
            var[:] = values
        for name, values in [('L11_fromN', links), ('L11_toN', targets)]:
            ds.createVariable(name, 'i4', ('links',), fill_value=-9999)[:] = values


def test_routing(tmp_path):
    for tool in ('gfortran', 'nf-config', 'make'):
        if not shutil.which(tool):
            pytest.skip(f'{tool} required')
    build = tmp_path / 'build'
    nf = shlex.split(subprocess.check_output(['nf-config', '--fflags'], text=True))
    flags = '-fimplicit-none -fcoarray=single -fdefault-real-8 -O0 -g -fcheck=all -J' + str(build / 'mods') + ' ' + ' '.join(nf)
    result = subprocess.run(['make', f'BUILD={build}', f'FFLAGS={flags}', str(build / 'objs/ncrouting.o'),
                             str(build / 'objs/netcdfflux.o')], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    libs = shlex.split(subprocess.check_output(['nf-config', '--flibs'], text=True))
    libs.insert(0, '-L' + subprocess.check_output(['nc-config', '--libdir'], text=True).strip())
    executable = tmp_path / 'test-routing'
    objects = list((build / 'objs').glob('*.o'))
    result = subprocess.run(['gfortran', '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
                             '-I', str(build / 'mods'), *nf,
                             str(ROOT / 'tests/fortran/test_ncrouting.f90'),
                             *map(str, objects), *libs, '-o', str(executable)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    valid = tmp_path / 'valid.nc'
    fixture(valid)
    result = subprocess.run([str(executable), str(valid), 'synthetic'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Synthetic routing checks passed' in result.stdout
    errors = {
        'cycle': 'cycle in routing graph', 'duplicate-link': 'duplicate downstream link',
        'bad-target': 'missing node', 'missing-link': 'no downstream link',
        'self-loop': 'self-loop', 'outlet-link': 'no downstream link',
        'duplicate-id': 'duplicate L11 ID', 'invalid-id': 'IDs must cover',
        'invalid-mask': 'invalid L11_domain_mask', 'invalid-code': 'invalid L11_fDir',
        'grid-mismatch': 'do not match Qrouted grid', 'nan': 'nonfinite',
        'missing-variable': 'NetCDF error', 'swapped-variable-dims': 'NetCDF error',
    }
    for label, expected in errors.items():
        path = tmp_path / (label + '.nc')
        fixture(path, label)
        result = subprocess.run([str(executable), str(path), 'invalid'], cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert expected in result.stdout, result.stdout + result.stderr
    result = subprocess.run([str(executable), 'unused', 'required'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'ADEnc requires drutes.conf/netcdf/mRM_restart_001.nc' in result.stderr + result.stdout
    inputs = tmp_path / 'drutes.conf/netcdf'
    inputs.mkdir(parents=True)
    (inputs / 'mRM_restart_001.nc').symlink_to(valid)
    result = subprocess.run([str(executable), 'unused', 'required'], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Mandatory routing read passed' in result.stdout
    # Optional local real-input regression; CI may not distribute large mHM data.
    routing = ROOT / 'drutes.conf/netcdf/mRM_restart_001.nc'
    forcing = ROOT / 'drutes.conf/netcdf/mRM_Fluxes_States.nc'
    if routing.exists() and forcing.exists():
        result = subprocess.run([str(executable), str(routing), 'real', str(forcing)],
                                cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        rows = np.genfromtxt(tmp_path / 'routing-analysis.csv', delimiter=',', names=True)
        nc = pytest.importorskip('netCDF4')
        with nc.Dataset(routing) as network, nc.Dataset(forcing) as data:
            ids = network['L11_Id'][:]
            mask = network['L11_domain_mask'][:]
            lat, lon = network['L11_domain_lat'][:], network['L11_domain_lon'][:]
            mapping = {}
            raw_codes = {}
            for row, col in zip(*np.where(mask == 1)):
                iy = int(np.argmin(abs(data['lat'][:] - lat[row, col])))
                ix = int(np.argmin(abs(data['lon'][:] - lon[row, col])))
                cell = iy * len(data['lon']) + ix + 1
                mapping[int(ids[row, col])] = cell
                raw_codes[cell] = int(network['L11_fDir'][row, col])
            expected = {cell: 0 for cell in mapping.values()}
            for source, target in zip(network['L11_fromN'][:].compressed(), network['L11_toN'][:].compressed()):
                expected[mapping[int(source)]] = mapping[int(target)]
            upstream = {cell: 0 for cell in expected}
            for target in expected.values():
                if target:
                    upstream[target] += 1
            assert len(rows) == len(mapping)
            for row in rows:
                cell = int(row['cell'])
                assert row['downstream_cell'] == expected[cell]
                assert row['upstream_count'] == upstream[cell]
                assert row['raw_fdir'] == raw_codes[cell]
                assert np.isclose(np.hypot(row['unit_x'], row['unit_y']), 1 if expected[cell] else 0)
