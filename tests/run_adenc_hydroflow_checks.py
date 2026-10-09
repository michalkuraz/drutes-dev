#!/usr/bin/env python3
"""Linux/Mac standalone acceptance checks; no pytest and never DRUtES main.

Build project objects first with bounds checking. All fixtures/results are
created in a NEW directory, never in existing simulation outputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output-parent", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = args.build.resolve()
    args.output_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="hydro-checks-", dir=args.output_parent.resolve()))
    objects = sorted(p for p in (build / "objs").glob("*.o") if p.name != "main.o")
    if not objects or not (build / "mods/nchydroflow.mod").exists():
        raise SystemExit("Build the current project with NetCDF before running these checks")
    nf_flags = shlex.split(subprocess.check_output(["nf-config", "--fflags"], text=True))
    libraries = shlex.split(subprocess.check_output(["nf-config", "--flibs"], text=True))
    libraries.insert(0, "-L" + subprocess.check_output(["nc-config", "--libdir"], text=True).strip())
    results: list[dict[str, object]] = []

    def run(label: str, command: list[str], cwd: Path, expected: str, fail: bool = False) -> None:
        result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=240)
        log = work / (label + ".log")
        log.write_text(result.stdout + result.stderr)
        passed = (result.returncode != 0 if fail else result.returncode == 0) and expected in log.read_text()
        results.append({"check": label, "passed": passed, "exit": result.returncode, "log": str(log)})
        (work / "RESULTS.json").write_text(json.dumps(results, indent=2) + "\n")
        if not passed:
            raise RuntimeError(f"{label} failed; inspect {log}\n{log.read_text()[-6000:]}")
        print(f"PASS {label}", flush=True)

    executables = {}
    for model in ("hydroflow", "conservative", "banks"):
        binary = work / ("test-" + model)
        command = ["gfortran", "-fcoarray=single", "-fdefault-real-8", "-fcheck=all",
                   "-I", str(build / "mods"), *nf_flags,
                   str(root / f"tests/fortran/test_adenc_{model}.f90"),
                   *map(str, objects), *libraries, "-o", str(binary)]
        run("compile-" + model, command, work, "")
        executables[model] = str(binary)
    for model, expected in (("conservative", "ADEnc conservative checks passed"),
                            ("banks", "closed mass checks passed")):
        case = work / model
        (case / "out").mkdir(parents=True)
        run("legacy-" + model, [executables[model]], case, expected)
    run("kernel", [executables["hydroflow"]], work, "Hydroflow kernel checks passed")
    config = "y\n0\n0\n1e-12\n1000\n100\n4\n1 3 -1 0\n3 5 -1 0\n2 4 1 0\n4 6 1 0\n"
    cases = {
        "daily": (config, "production", "Hydroflow production checks passed", False),
        "linear": (config.replace("y\n0\n0\n", "y\n0\n1\n", 1), "linear",
                   "Hydroflow production checks passed", False),
        "off": ("n\n", "off", "Hydroflow config checks passed", False),
        "absent": (None, "absent", "Hydroflow config checks passed", False),
        "duplicate": (config.replace("3 5 -1 0", "1 3 -1 0"), "production", "Duplicate hydroflow port", True),
        "interior": (config.replace("1 3 -1 0", "1 4 -1 0"), "production", "active-domain boundary", True),
        "correction": (config.replace("\n100\n4\n", "\n0.001\n4\n"), "production", "correction exceeds", True),
        "source": (config.replace("y\n0\n", "y\n2\n", 1), "production", "policy must be 0 or 1", True),
        "legacy-mode": (config, "legacy", "requires conservative", True),
    }
    for label, (content, argument, expected, fail) in cases.items():
        case = work / label
        (case / "out").mkdir(parents=True)
        inputs = case / "drutes.conf/netcdf"
        inputs.mkdir(parents=True)
        if content is not None:
            (inputs / "hydroflow.conf").write_text(content)
        run(label, [executables["hydroflow"], argument], case, expected, fail)
    lateral_config = config.replace("y\n0\n0\n", "y\n1\n1\n", 1)
    # Two disjoint groups, each TOTAL Q=.1 initially, increasing to .2.
    # Physical volumes/load are distributed per area, never duplicated per FE.
    def lateral_data(concentration: int = 1, zero: bool = False) -> str:
        q0, q1 = (0, 0) if zero else (.1, .2)
        series = (f"0 {q0} {concentration}\n3600 {q1} {concentration}\n"
                  f"172800 {q1} {concentration}\n259200 {q1} {concentration}\n")
        return "2\n2 4\n1\n2\n" + series + "2 4\n3\n4\n" + series

    lateral_cases = {
        "lateral-uniform": (lateral_data(), "lateral-uniform", "production checks passed", False),
        "lateral-clean": (lateral_data(0), "lateral-clean", "production checks passed", False),
        "lateral-zero": (lateral_data(1, True), "lateral-zero", "production checks passed", False),
        "lateral-load": (lateral_data().replace("3600 0.2 1", "3600 0.2 2"), "lateral-load", "production checks passed", False),
        "lateral-overlap": (lateral_data().replace("3\n4\n", "1\n2\n"), "lateral-overlap", "production checks passed", False),
        "lateral-missing": (None, "production", "policy 1 requires lateral.conf", True),
        "lateral-negative": (lateral_data().replace("0 0.1 1", "0 -0.1 1"), "production", "nonnegative", True),
        "lateral-nan": (lateral_data().replace("0 0.1 1", "0 NaN 1"), "production", "Nonfinite lateral", True),
        "lateral-duplicate": (lateral_data().replace("1\n2\n", "1\n1\n", 1), "production", "Duplicate lateral", True),
        "lateral-element": (lateral_data().replace("1\n2\n", "1\n99\n", 1), "production", "outside FE mesh", True),
        "lateral-time": (lateral_data().replace("3600", "0", 1), "production", "strictly increase", True),
        "lateral-short": (lateral_data().replace("259200", "172801"), "lateral-short", "does not cover", True),
        "lateral-coverage": (lateral_data().replace("172800", "7200").replace("259200", "10800"), "production", "does not cover configured", True),
        "lateral-extra": (lateral_data()+"99\n", "production", "Unexpected trailing lateral", True),
        "lateral-extra-field": (lateral_data().replace("0 0.1 1", "0 0.1 1 9"), "production", "Extra or invalid", True),
    }
    for label, (content, argument, expected, fail) in lateral_cases.items():
        case = work / label
        inputs = case / 'drutes.conf/netcdf'
        inputs.mkdir(parents=True)
        (case / 'out').mkdir()
        (inputs / 'hydroflow.conf').write_text(lateral_config)
        if content is not None:
            (inputs / 'lateral.conf').write_text(content)
        run(label, [executables['hydroflow'], argument], case, expected, fail)
    print(f"All {len(results)} checks passed. Preserved logs: {work}")


if __name__ == "__main__":
    main()
