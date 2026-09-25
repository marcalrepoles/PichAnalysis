"""Resolve application-owned resources independently of the process CWD."""
from __future__ import annotations

import os
import sys
from pathlib import Path


def r_scripts_dir(*, required: bool = True) -> Path:
    candidates: list[Path] = []
    override = os.environ.get("PICHANALYSIS_RESOURCE_ROOT")
    if override:
        candidates.append(Path(override) / "r_scripts")
    frozen_root = getattr(sys, "_MEIPASS", None)
    if frozen_root:
        candidates.append(Path(frozen_root) / "r_scripts")
    if getattr(sys, "frozen", False):
        candidates.append(Path(sys.executable).resolve().parent / "r_scripts")
    if not getattr(sys, "frozen", False):
        candidates.append(Path(__file__).resolve().parents[3] / "r_scripts")
    for folder in candidates:
        if (folder / "00_runtime_test.R").is_file():
            return folder.resolve()
    if not required:
        return candidates[0].resolve()
    raise FileNotFoundError(
        "PichAnalysis R scripts were not found. Check the application installation "
        "or PICHANALYSIS_RESOURCE_ROOT."
    )


def r_script(name: str) -> Path:
    if Path(name).name != name or not name.endswith(".R"):
        raise ValueError("Invalid R script name.")
    path = r_scripts_dir() / name
    if not path.is_file():
        raise FileNotFoundError(f"Required PichAnalysis R script is missing: {name}")
    return path

def branding_asset(name: str) -> Path:
    """Locate an application-owned image in source or frozen builds."""
    if name not in {"icon.ico", "splash.png"}:
        raise ValueError("Invalid branding asset name.")
    frozen_root = getattr(sys, "_MEIPASS", None)
    root = Path(frozen_root) if frozen_root else Path(__file__).resolve().parents[3]
    path = root / "branding" / name if frozen_root else root / name
    if not path.is_file():
        raise FileNotFoundError(f"Required PichAnalysis branding asset is missing: {name}")
    return path
