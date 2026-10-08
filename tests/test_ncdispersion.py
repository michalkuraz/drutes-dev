"""Isolated checks of production dispersion parsing and tensor algebra; no model run."""
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_dispersion(tmp_path: Path) -> None:
    compiler = shutil.which('gfortran')
    if compiler is None:
        pytest.skip('gfortran required')
    executable = tmp_path / 'test-dispersion'
    result = subprocess.run([
        compiler, '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-O0',
        '-J', str(tmp_path), '-I', str(tmp_path),
        str(ROOT / 'src/core/typy.f90'),
        str(ROOT / 'src/models/fluxLS/ncdispersion.f90'),
        str(ROOT / 'tests/fortran/test_ncdispersion.f90'), '-o', str(executable),
    ], cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([str(executable)], cwd=tmp_path, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'ADEnc dispersion checks passed' in result.stdout
