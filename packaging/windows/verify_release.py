"""Inspect and execute a real Windows onedir or installed bundle."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path


CRITICAL_R = tuple(f"{index:02d}_{name}.R" for index, name in enumerate((
    "runtime_test", "mapping_annotation", "presence_absence", "go_analysis",
    "kegg_analysis", "reactome_analysis", "mitocarta_analysis",
    "interpro_pfam_analysis", "string_analysis", "complex_analysis",
    "mtdna_analysis", "proteomics_qc", "differential_preparation",
    "differential_analysis",
)))


def verify(bundle: Path, output: Path, *, missing_r: bool = True) -> dict:
    bundle, output = bundle.resolve(), output.resolve()
    exe = bundle / "PichAnalysis.exe"
    if not exe.is_file():
        raise AssertionError(f"Missing Windows executable: {exe}")
    internal = bundle / "_internal"
    scripts = internal / "r_scripts"
    if not scripts.is_dir():
        scripts = bundle / "r_scripts"
    for name in CRITICAL_R:
        if not (scripts / name).is_file():
            raise AssertionError(f"Missing bundled R script: {name}")
    if not (scripts / "lib" / "common.R").is_file():
        raise AssertionError("Missing bundled R helpers")
    for forbidden in (bundle / "tests", internal / "tests", scripts / "tests", bundle / ".git"):
        if forbidden.exists():
            raise AssertionError(f"Forbidden distribution content: {forbidden}")
    for forbidden in ("GO", "KEGG", "Reactome", "MitoCarta", "STRING", "ComplexPortal", "InterPro"):
        if (bundle / "databases" / forbidden).exists():
            raise AssertionError(f"Scientific database snapshot bundled: {forbidden}")
    plugins = list(bundle.rglob("qwindows.dll"))
    if not plugins:
        raise AssertionError("Qt qwindows platform plugin is missing")
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    from pichanalysis.core.r_runtime import RRuntime
    external_r = RRuntime()
    if external_r.available:
        library = subprocess.run((external_r.executable, "--vanilla", "-e", "cat(.libPaths()[1])"),
            capture_output=True, text=True, timeout=15, check=False)
        if library.returncode == 0 and library.stdout.strip():
            environment["R_LIBS_USER"] = library.stdout.strip()
    environment["QT_QPA_PLATFORM"] = "offscreen"
    environment.pop("PICHANALYSIS_RESOURCE_ROOT", None)
    environment["LOCALAPPDATA"] = str(output / "profile" / "Local")
    environment["APPDATA"] = str(output / "profile" / "Roaming")
    environment["PICHANALYSIS_BOOT_ERROR_FILE"] = str(output / "bootstrap_error.txt")
    system_root = Path(environment.get("SystemRoot", r"C:\Windows"))
    environment["PATH"] = os.pathsep.join((str(system_root / "System32"), str(system_root)))
    cwd = Path(tempfile.mkdtemp(prefix="pich_release_cwd_"))
    checks = []
    try:
        variants = (False, True) if missing_r else (False,)
        for without_r in variants:
            command = [str(exe), "--release-smoke", str(output)]
            if without_r:
                command.append("--simulate-missing-r")
            process = subprocess.run(command, cwd=cwd, env=environment, timeout=120, check=False)
            report = output / ("smoke_missing_r.json" if without_r else "smoke.json")
            if not report.is_file():
                raise AssertionError(f"Smoke produced no result file (exit {process.returncode}): {report}")
            detail = json.loads(report.read_text(encoding="utf-8"))
            if process.returncode or detail.get("status") != "passed":
                raise AssertionError(f"Packaged smoke failed: {detail}")
            if not detail["preflight"]["resources"] or not detail["preflight"]["python"]:
                raise AssertionError(f"Packaged resource/import preflight failed: {detail}")
            if without_r and detail["preflight"]["rscript"]:
                raise AssertionError("Missing-R simulation failed")
            if not str(detail["r_scripts_dir"]).lower().startswith(str(bundle).lower()):
                raise AssertionError("Frozen resource lookup escaped the inspected bundle")
            if not without_r and not detail["preflight"]["required_r_packages"]:
                raise AssertionError(f"R package preflight failed: {detail}")
            if not without_r and not detail.get("r_resource_execution"):
                raise AssertionError("Bundled R script was not executed")
            checks.append(detail)
    finally:
        cwd.rmdir()
    return {"bundle": str(bundle), "qt_plugins": [str(path.relative_to(bundle)) for path in plugins],
            "checks": checks, "status": "passed"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.bundle, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
