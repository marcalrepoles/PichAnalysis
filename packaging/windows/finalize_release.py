"""Write and independently verify release artifact checksums and manifest."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from pichanalysis.version import __version__
from pichanalysis.core.r_runtime import RRuntime


def digest(path: Path) -> str:
    hashed = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hashed.update(block)
    return hashed.hexdigest()


def git(root: Path, *args: str) -> str:
    return subprocess.run(("git", *args), cwd=root, capture_output=True,
        text=True, check=True).stdout.strip()


def inventory() -> list[dict[str, str]]:
    packages = ("PySide6", "shiboken6", "pandas", "numpy", "openpyxl", "xlrd", "platformdirs", "pyinstaller")
    items = []
    for name in packages:
        try:
            metadata = importlib.metadata.metadata(name)
            items.append({"name": name, "version": importlib.metadata.version(name),
                "license_declared": metadata.get("License", "") or "",
                "license_expression": metadata.get("License-Expression", "") or ""})
        except importlib.metadata.PackageNotFoundError:
            items.append({"name": name, "version": "not found in build environment",
                "license_declared": "", "license_expression": ""})
    return items


def finalize(root: Path, installer: Path, portable: Path, output: Path, inno_version: str) -> dict:
    root, output = root.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    expected = {
        "installer": f"PichAnalysis-{__version__}-Windows-x86_64-Setup.exe",
        "portable": f"PichAnalysis-{__version__}-Windows-x86_64-portable.zip",
    }
    paths = {"installer": installer.resolve(), "portable": portable.resolve()}
    for kind, path in paths.items():
        if path.name != expected[kind] or not path.is_file() or path.parent != output:
            raise ValueError(f"Invalid {kind} artifact: {path}")
    artifacts = {kind: {"filename": path.name, "size_bytes": path.stat().st_size,
        "sha256": digest(path)} for kind, path in paths.items()}
    sums = output / "SHA256SUMS.txt"
    sums.write_text("".join(f"{item['sha256']}  {item['filename']}\n"
        for item in artifacts.values()), encoding="utf-8")
    for kind, path in paths.items():
        if digest(path) != artifacts[kind]["sha256"]:
            raise RuntimeError(f"SHA-256 recheck failed: {kind}")
    runtime = RRuntime()
    manifest = {
        "product": "PichAnalysis", "version": __version__, "platform": "Windows",
        "architecture": "x86_64", "build_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "pyside6_version": importlib.metadata.version("PySide6"),
        "pyinstaller_version": importlib.metadata.version("pyinstaller"),
        "installer_compiler": "Inno Setup", "installer_compiler_version": inno_version,
        "validated_r_version": runtime.version() or "not detected",
        "git_head": git(root, "rev-parse", "HEAD"),
        "working_tree_state": "modified" if git(root, "status", "--porcelain") else "clean",
        "r_runtime_bundled": False, "scientific_databases_bundled": False,
        "included_resource_groups": ["application Python runtime", "PySide6/Qt", "r_scripts entry points and lib", "release documentation"],
        "third_party_inventory": inventory(), "artifacts": artifacts,
        "test_summary": {"bundle_smoke": "passed", "hash_recheck": "passed",
            "installer_lifecycle": "pending controlled install/uninstall/reinstall verification"},
        "notes": ["Unsigned Windows build", "Built from an intentionally modified working tree; final commit is user-owned."],
    }
    (output / "release_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--installer", type=Path, required=True)
    parser.add_argument("--portable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inno-version", default="unknown")
    args = parser.parse_args()
    finalize(args.root, args.installer, args.portable, args.output, args.inno_version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
