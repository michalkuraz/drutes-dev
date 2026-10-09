"""Compile and exercise the production Fortran width routines in isolated directories."""

from pathlib import Path
import shlex
import shutil
import subprocess
import numpy as np

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], cwd: Path) -> str:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=240)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def test_contact_geometry(tmp_path: Path) -> None:
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for Fortran geometry tests")
    executable = tmp_path / "test-width"
    run(
        [compiler, "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-O0",
         "-J", str(tmp_path), "-I", str(tmp_path),
         str(ROOT / "src/core/typy.f90"),
         str(ROOT / "src/models/fluxLS/ncwidth_geometry.f90"),
         str(ROOT / "tests/fortran/test_ncwidth_geometry.f90"),
         "-o", str(executable)],
        tmp_path,
    )
    assert "river contact geometry checks" in run([str(executable)], tmp_path)


def test_width_cache_with_real_modules(tmp_path: Path) -> None:
    required = {name: shutil.which(name) for name in ("gfortran", "make", "nf-config", "nc-config")}
    if not all(required.values()):
        pytest.skip("Fortran, Make and NetCDF development tools are required")
    build = tmp_path / "build"
    binary_directory = tmp_path / "bin"
    # Build only; never run bin/drutes (its main program deletes out/*).
    # Bounds-check the changed modules while keeping the real project interfaces.
    nf_flags = shlex.split(run([required["nf-config"], "--fflags"], tmp_path))
    flags = " ".join([
        "-DHAVE_NETCDF -fimplicit-none -fcoarray=single -fdefault-real-8",
        "-O0 -g -fcheck=all -J" + shlex.quote(str(build / "mods")),
        *map(shlex.quote, nf_flags),
    ])
    run([required["make"], "build_target", f"BUILD={build}", f"BINDIR={binary_directory}",
         f"FFLAGS={flags}"], ROOT)
    objects = sorted(path for path in (build / "objs").glob("*.o") if path.name != "main.o")
    assert objects
    executable = tmp_path / "test-width-cache"
    nc_libdir = run([required["nc-config"], "--libdir"], tmp_path).strip()
    libraries = [f"-L{nc_libdir}", *shlex.split(run([required["nf-config"], "--flibs"], tmp_path))]
    run([required["gfortran"], "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
         "-I", str(build / "mods"), *nf_flags,
         str(ROOT / "tests/fortran/test_ncflux_widths.f90"),
         *map(str, objects), *libraries, "-o", str(executable)], tmp_path)
    assert "river width integration checks" in run([str(executable)], tmp_path)

    dispersion_executable = tmp_path / "test-dispersion-callback"
    run([required["gfortran"], "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
         "-I", str(build / "mods"), *nf_flags,
         str(ROOT / "tests/fortran/test_adenc_dispersion_callback.f90"),
         *map(str, objects), *libraries, "-o", str(dispersion_executable)], tmp_path)
    assert "ADEnc dispersion callback checks passed" in run([str(dispersion_executable)], tmp_path)

    supg_executable = tmp_path / "test-adenc-supg"
    run([required["gfortran"], "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
         "-I", str(build / "mods"), *nf_flags,
         str(ROOT / "tests/fortran/test_adenc_supg.f90"),
         *map(str, objects), *libraries, "-o", str(supg_executable)], tmp_path)
    assert "ADEnc SUPG checks passed" in run([str(supg_executable)], tmp_path)
    bank_executable = tmp_path / 'test-adenc-banks'
    run([required['gfortran'], '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
         '-I', str(build / 'mods'), *nf_flags,
         str(ROOT / 'tests/fortran/test_adenc_banks.f90'),
         *map(str, objects), *libraries, '-o', str(bank_executable)], tmp_path)
    assert 'closed mass checks passed' in run([str(bank_executable)], tmp_path)
    conservative_executable = tmp_path / 'test-adenc-conservative'
    run([required['gfortran'], '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
         '-I', str(build / 'mods'), *nf_flags,
         str(ROOT / 'tests/fortran/test_adenc_conservative.f90'),
         *map(str, objects), *libraries, '-o', str(conservative_executable)], tmp_path)
    (tmp_path / 'out').mkdir(exist_ok=True)
    assert 'ADEnc conservative checks passed' in run([str(conservative_executable)], tmp_path)
    budget = np.genfromtxt(tmp_path / 'out/adenc_mass_balance.csv', delimiter=',', names=True)
    assert budget.size == 20  # accepted pulse/outflow steps; rejected trials add no rows
    assert np.all(np.isfinite(budget.tolist()))
    assert np.max(np.abs(budget['error'])) < 1e-11
    assert np.max(np.abs(budget['free_residual'])) < 1e-11
    assert budget['cumulative_in'][-1] > 0
    assert budget['cumulative_out'][-1] > 0
    assert budget['cumulative_source'][-1] == 0
    hydro_executable = tmp_path / 'test-adenc-hydroflow'
    run([required['gfortran'], '-fcoarray=single', '-fdefault-real-8', '-fcheck=all',
         '-I', str(build / 'mods'), *nf_flags,
         str(ROOT / 'tests/fortran/test_adenc_hydroflow.f90'),
         *map(str, objects), *libraries, '-o', str(hydro_executable)], tmp_path)
    assert 'Hydroflow kernel checks passed' in run([str(hydro_executable)], tmp_path)
    hydro_config = 'y\n0\n0\n1e-12\n1000\n100\n4\n1 3 -1 0\n3 5 -1 0\n2 4 1 0\n4 6 1 0\n'
    hydro_case = tmp_path / 'hydro-production'
    hydro_inputs = hydro_case / 'drutes.conf/netcdf'
    hydro_inputs.mkdir(parents=True)
    (hydro_case / 'out').mkdir()
    (hydro_inputs / 'hydroflow.conf').write_text(hydro_config)
    assert 'Hydroflow production checks passed' in run([str(hydro_executable), 'production'], hydro_case)
    linear_case = tmp_path / 'hydro-linear'
    linear_inputs = linear_case / 'drutes.conf/netcdf'
    linear_inputs.mkdir(parents=True)
    (linear_case / 'out').mkdir()
    (linear_inputs / 'hydroflow.conf').write_text(hydro_config.replace('y\n0\n0\n', 'y\n0\n1\n', 1))
    assert 'Hydroflow production checks passed' in run([str(hydro_executable), 'linear'], linear_case)
    for label, suffix, argument, succeeds in (
        ('enabled', 'y\n', 'island', True),
        ('point-contact', 'y\n', 'island-point', True),
        ('disabled', 'n\n', 'island', False),
        ('default', '', 'island', False),
    ):
        case = tmp_path / f'hydro-filter-{label}'
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        (case / 'out').mkdir()
        (inputs / 'hydroflow.conf').write_text(hydro_config + suffix)
        result = subprocess.run([str(hydro_executable), argument], cwd=case,
                                text=True, capture_output=True, timeout=20)
        output = result.stdout + result.stderr
        assert (result.returncode == 0) == succeeds, output
        assert ('component filter checks passed' if succeeds else
                'Every hydroflow component needs') in output, output
    invalid_hydro = {
        'legacy': (hydro_config, 'requires conservative'),
        'source': (hydro_config.replace('y\n0\n', 'y\n2\n', 1), 'policy must be 0 or 1'),
        'nan': (hydro_config.replace('1e-12', 'NaN'), 'Nonfinite hydroflow'),
        'duplicate': (hydro_config.replace('3 5 -1 0', '1 3 -1 0'), 'Duplicate hydroflow port'),
        'interior': (hydro_config.replace('1 3 -1 0', '1 4 -1 0'), 'active-domain boundary'),
        'no-outlet': (hydro_config.replace('2 4 1 0', '2 4 -1 0').replace('4 6 1 0', '4 6 -1 0'),
                      'real concentration boundary ID'),
        'correction': (hydro_config.replace('\n100\n4\n', '\n0.001\n4\n'), 'correction exceeds'),
        'extra': (hydro_config.replace('1 3 -1 0', '1 3 -1 0 junk'), 'Extra or invalid'),
    }
    for label, (content, expected) in invalid_hydro.items():
        case = tmp_path / f'hydro-{label}'
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        (case / 'out').mkdir()
        (inputs / 'hydroflow.conf').write_text(content)
        result = subprocess.run([str(hydro_executable), 'legacy' if label == 'legacy' else 'production'],
                                cwd=case, text=True, capture_output=True, timeout=20)
        assert result.returncode != 0, label
        assert expected in result.stdout + result.stderr, result.stdout + result.stderr
    for label, content in {'absent': None, 'off': 'n\n'}.items():
        case = tmp_path / f'hydro-{label}'
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        if content is not None:
            (inputs / 'hydroflow.conf').write_text(content)
        assert 'Hydroflow config checks passed' in run([str(hydro_executable), label], case)
    for label, content in {'absent': None, 'off': 'n\nn\n', 'on': 'y\ny\n',
                           'audit': 'n\ny\n', 'schwarz': 'y\ny\n', '3d': 'y\ny\n',
                           'backup': 'y\ny\n'}.items():
        case = tmp_path / f'conservative-{label}'
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        if content is not None:
            (inputs / 'conservative.conf').write_text(content)
        result = subprocess.run([str(conservative_executable), label], cwd=case,
                                text=True, capture_output=True, timeout=20)
        if label in ('schwarz', '3d', 'backup'):
            assert result.returncode != 0
            assert ('requires 2D standard Picard' if label != 'backup' else
                    'restart from backup is not supported') in result.stdout + result.stderr
        else:
            assert result.returncode == 0, result.stdout + result.stderr
    case = tmp_path / 'conservative-inflow'
    case.mkdir()
    (case / 'out').mkdir()
    result = subprocess.run([str(conservative_executable), 'inflow'], cwd=case,
                            text=True, capture_output=True, timeout=20)
    assert result.returncode != 0
    assert 'Label exterior inflow' in result.stdout + result.stderr
    for label, content in {'absent': None, 'on': 'y\n', 'schwarz': 'y\n', '3d': 'y\n'}.items():
        case = tmp_path / f'banks-{label}'
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        if content is not None:
            (inputs / 'riverbank.conf').write_text(content)
        result = subprocess.run([str(bank_executable), label], cwd=case,
                                text=True, capture_output=True, timeout=20)
        if label in ('schwarz', '3d'):
            assert result.returncode != 0
            assert 'require 2D and standard Picard' in result.stdout + result.stderr
        else:
            assert result.returncode == 0, result.stdout + result.stderr
    for label, content in {"absent": None, "off": "n\n1.0\n", "on": "y\n1.0\n"}.items():
        case = tmp_path / f"supg-{label}"
        inputs = case / "drutes.conf/netcdf"
        inputs.mkdir(parents=True)
        (case / "out").mkdir()
        if content is not None:
            (inputs / "supg.conf").write_text(content)
        assert "SUPG configuration checks passed" in run([str(supg_executable), label], case)
    for label, factor in {"negative": "-1", "nan": "NaN", "infinite": "Inf"}.items():
        case = tmp_path / f"supg-{label}"
        inputs = case / "drutes.conf/netcdf"
        inputs.mkdir(parents=True)
        (case / "out").mkdir()
        (inputs / "supg.conf").write_text(f"y\n{factor}\n")
        result = subprocess.run([str(supg_executable), "invalid"], cwd=case,
                                text=True, capture_output=True, timeout=20)
        assert result.returncode != 0
        assert "SUPG factor must be" in result.stdout + result.stderr

    for label, content in {"shock-absent": None, "shock-off": "n\n1\n", "shock": "y\n1\n"}.items():
        case = tmp_path / label
        inputs = case / "drutes.conf/netcdf"
        inputs.mkdir(parents=True)
        (case / "out").mkdir()
        if content is not None:
            (inputs / "shock.conf").write_text(content)
        assert "SUPG configuration checks passed" in run([str(supg_executable), label], case)
    for label, factor in {"negative": "-1", "nan": "NaN", "infinite": "Inf"}.items():
        case = tmp_path / f"shock-invalid-{label}"
        inputs = case / "drutes.conf/netcdf"
        inputs.mkdir(parents=True)
        (case / "out").mkdir()
        (inputs / "shock.conf").write_text(f"y\n{factor}\n")
        result = subprocess.run([str(supg_executable), "invalid"], cwd=case,
                                text=True, capture_output=True, timeout=20)
        assert result.returncode != 0
        assert "Shock factor must be" in result.stdout + result.stderr

    # Real local assembly of a tiny synthetic convection-dominated strip.
    # Does not call DRUtES main or use any repository simulation inputs/outputs.
    solutions = {}
    for mode in ("strip-off", "strip-on"):
        case = tmp_path / mode
        case.mkdir()
        assert "SUPG strip assembled" in run([str(supg_executable), mode], case)
        matrix = np.loadtxt(case / "benchmark-matrix.dat")
        rhs = np.loadtxt(case / "benchmark-rhs.dat")
        solutions[mode] = np.linalg.solve(matrix, rhs)
        assert np.max(np.abs(matrix @ solutions[mode] - rhs)) < 1e-10
    assert solutions["strip-off"].min() < -0.1  # Demonstrate unstabilized oscillations.
    # SUPG is not a positivity-preserving limiter on this anisotropic triangle mesh.
    # Test its actual contract: strongly reduced oscillation and solution error.
    assert abs(solutions["strip-on"].min()) < 0.2 * abs(solutions["strip-off"].min())
    assert solutions["strip-on"].max() <= 1 + 1e-10
    x = np.tile(np.linspace(0, 1, 21), 2)
    exact = (np.exp((x - 1) / 0.005) - np.exp(-1 / 0.005)) / (1 - np.exp(-1 / 0.005))
    rms = {mode: np.sqrt(np.mean((values - exact)**2)) for mode, values in solutions.items()}
    assert rms["strip-on"] < 0.3 * rms["strip-off"]

    # Nonlinear Picard solves: viscosity uses the preceding iterate, not a
    # manufactured exact field. Underrelax this deliberately skinny strip.
    iterate = solutions["strip-on"].copy()
    for iteration in range(100):
        case = tmp_path / f"strip-shock-{iteration}"
        case.mkdir()
        np.savetxt(case / "iterate.dat", iterate[None, :])
        run([str(supg_executable), "strip-shock"], case)
        matrix = np.loadtxt(case / "benchmark-matrix.dat")
        rhs = np.loadtxt(case / "benchmark-rhs.dat")
        solution = np.linalg.solve(matrix, rhs)
        error = np.max(np.abs(solution - iterate))
        if error < 1e-9:
            break
        iterate = 0.5 * (iterate + solution)
    assert error < 1e-9, f"shock Picard did not converge: {error}"
    assert abs(solution.min()) < abs(solutions["strip-on"].min())
    assert solution.max() <= 1 + 1e-10
    # Additional viscosity reduces undershoot but can smear a boundary layer;
    # do not demand or claim better exact error than SUPG alone.
    assert np.sqrt(np.mean((solution-exact)**2)) < rms["strip-off"]
    print("Shock strip:", iteration + 1, solution.min(), solution.max(),
          np.sqrt(np.mean((solution-exact)**2)))

    channel_executable = tmp_path / "test-channel-paths"
    run([required["gfortran"], "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
         "-I", str(build / "mods"), *nf_flags,
         str(ROOT / "tests/fortran/test_ncchannel_paths.f90"),
         *map(str, objects), *libraries, "-o", str(channel_executable)], tmp_path)

    for channel_count in (1, 2, 5):
        case_directory = tmp_path / f"channels-{channel_count}"
        input_directory = case_directory / "drutes.conf" / "netcdf"
        input_directory.mkdir(parents=True)
        (input_directory / "channel.dat").write_text(
            "# first channel\n0.0 0.0\n1.0 0.0\n2.0 0.0\n",
            encoding="utf-8",
        )
        for index in range(2, channel_count + 1):
            x = 100.0 * index
            (input_directory / f"channel{index}.dat").write_text(
                f"# channel {index}\n{x} {index}.0\n{x + 1.0} {index}.0\n",
                encoding="utf-8",
            )
        assert "multiple channel path checks" in run(
            [str(channel_executable), str(channel_count)], case_directory
        )

    # Bad records must fail instead of silently truncating the second channel.
    invalid = {
        "malformed": ("200 2\n201 2\nbad coordinates\n202 2\n", 3),
        "missing-coordinate": ("200 2\n201 2\n202\n203 2\n", 3),
        "extra-coordinate": ("200 2\n201 2\n202 2 9\n", 3),
        "nan": ("200 2\n201 2\nNaN 2\n", 3),
        "inf": ("200 2\n201 2\n202 Inf\n", 3),
        "duplicate": ("200 2\n200 2\n201 2\n", 2),
        "null-coordinate": ("200 2\n201 2\n, 2\n", 3),
    }
    for label, (content, line) in invalid.items():
        case_directory = tmp_path / label
        inputs = case_directory / "drutes.conf" / "netcdf"
        inputs.mkdir(parents=True)
        (inputs / "channel.dat").write_text("0 0\n1 0\n2 0\n", encoding="utf-8")
        (inputs / "channel2.dat").write_text(content, encoding="utf-8")
        result = subprocess.run([str(channel_executable), "2"], cwd=case_directory,
                                text=True, capture_output=True, timeout=20)
        assert result.returncode != 0, label
        assert f"channel2.dat at line {line}:" in result.stdout, result.stdout + result.stderr

    # Sharing a confluence point between files is legal, with no joining edge.
    inputs = tmp_path / "confluence" / "drutes.conf" / "netcdf"
    inputs.mkdir(parents=True)
    (inputs / "channel.dat").write_text("0 0\n1 0\n200 2\n", encoding="utf-8")
    (inputs / "channel2.dat").write_text(
        "# branch\n\n200\t2 # confluence\n201 2\n# final comment\n", encoding="utf-8"
    )
    assert "multiple channel path checks" in run([str(channel_executable), "2"], inputs.parents[1])

    boundary_executable = tmp_path / "test-adenc-boundaries"
    run([required["gfortran"], "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
         "-I", str(build / "mods"), *nf_flags,
         str(ROOT / "tests/fortran/test_adenc_boundaries.f90"),
         *map(str, objects), *libraries, "-o", str(boundary_executable)], tmp_path)
    for inlets in (1, 2):
        case_directory = tmp_path / f"inlets-{inlets}"
        inputs = case_directory / "drutes.conf" / "netcdf"
        inputs.mkdir(parents=True)
        records = []
        for i in range(inlets):
            boundary_id = 101 + i
            records.append(f"{boundary_id} 1 y 0.0\n")
            (inputs / f"{boundary_id}.bc").write_text(
                f"0.0 {i+1}.0\n20.0 0.0\n", encoding="utf-8"
            )
        records.append(f"{101+inlets} 1 n 0.0\n")
        (case_directory / "boundaries.conf").write_text("".join(records), encoding="utf-8")
        assert "ADEnc boundary checks" in run([str(boundary_executable), str(inlets)], case_directory)

    bad_boundaries = tmp_path / "inlets-2" / "boundaries.conf"
    bad_boundaries.write_text("102 1 n 1.0\n101 1 n 2.0\n103 1 n 0.0\n", encoding="utf-8")
    result = subprocess.run([str(boundary_executable), "2"], cwd=bad_boundaries.parent,
                            text=True, capture_output=True, timeout=20)
    assert result.returncode != 0
    assert "expected ID" in result.stdout
