"""Read-only release readiness checks; never installs R or Python packages."""
from __future__ import annotations

import importlib.util
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .r_runtime import RRuntime
from .external_process import run_external
from .resources import r_scripts_dir


PYTHON_RUNTIME = ("PySide6", "pandas", "openpyxl", "xlrd", "platformdirs")
R_REQUIRED_BY = {
    "jsonlite": "Mapping, Presence/Absence, all downstream analyses",
    "readr": "Mapping, Presence/Absence, GO, KEGG, Reactome, MitoCarta, InterPro/Pfam",
    "openxlsx": "Mapping and scientific workbooks",
    "ggplot2": "Presence/Absence and scientific plots",
    "httr2": "network-backed identifier mapping",
    "DBI": "Mapping",
    "RSQLite": "Mapping",
    "AnnotationDbi": "GO",
    "GO.db": "GO",
    "xml2": "KEGG KGML",
    "igraph": "STRING",
    "limma": "Differential Statistics",
}
R_OPTIONAL_BY = {
    "org.Hs.eg.db": "GO, human projects only",
    "org.Mm.eg.db": "GO, mouse projects only",
    "MsCoreUtils": "Differential Preparation imputation",
    "imputeLCMD": "MinProb/QRILC imputation",
    "impute": "KNN imputation",
    "statmod": "robust limma moderation",
}
CRITICAL_SCRIPTS = tuple(f"{number:02d}_{name}.R" for number, name in enumerate((
    "runtime_test", "mapping_annotation", "presence_absence", "go_analysis",
    "kegg_analysis", "reactome_analysis", "mitocarta_analysis",
    "interpro_pfam_analysis", "string_analysis", "complex_analysis",
    "mtdna_analysis", "proteomics_qc", "differential_preparation",
    "differential_analysis",
)))


@dataclass(frozen=True)
class Check:
    name: str
    required_by: str
    installed: bool
    version: str = ""
    detail: str = ""


@dataclass(frozen=True)
class Preflight:
    python: tuple[Check, ...]
    rscript: Check
    r_packages: tuple[Check, ...]
    resources: tuple[Check, ...]
    writable: tuple[Check, ...]

    @property
    def ready(self) -> bool:
        return all(item.installed for item in (
            *self.python, self.rscript, *self.r_packages, *self.resources, *self.writable
        ) if item.name not in R_OPTIONAL_BY)


def _r_packages(executable: str) -> tuple[Check, ...]:
    names = tuple(R_REQUIRED_BY) + tuple(R_OPTIONAL_BY)
    quoted = ",".join('"' + name.replace('"', '\\"') + '"' for name in names)
    expression = (
        f"for(p in c({quoted})){{"
        "ok<-requireNamespace(p,quietly=TRUE);"
        "v<-if(ok)as.character(utils::packageVersion(p))else '';"
        "cat(p,'\\t',ok,'\\t',v,'\\n',sep='')}"
    )
    try:
        result = run_external((executable, "--vanilla", "-e", expression),
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=45, shell=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        return tuple(Check(name, R_REQUIRED_BY.get(name, R_OPTIONAL_BY.get(name, "")),
            False, detail=str(error)) for name in names)
    found = {}
    if result.returncode == 0:
        for line in result.stdout.splitlines():
            columns = line.split("\t")
            if len(columns) >= 3:
                found[columns[0]] = (columns[1] == "TRUE", columns[2])
    return tuple(Check(name, R_REQUIRED_BY.get(name, R_OPTIONAL_BY.get(name, "")),
        found.get(name, (False, ""))[0], found.get(name, (False, ""))[1],
        "" if result.returncode == 0 else (result.stderr or "R package check failed"))
        for name in names)


def _writable(directory: Path, label: str) -> Check:
    try:
        with tempfile.TemporaryFile(dir=directory):
            pass
    except OSError as error:
        return Check(label, "Application output", False, detail=str(error))
    return Check(label, "Application output", True)


def run_preflight(project_directory: Path | None = None,
                  runtime: RRuntime | None = None) -> Preflight:
    python = tuple(Check(name, "Python application runtime",
        importlib.util.find_spec(name) is not None) for name in PYTHON_RUNTIME)
    runtime = runtime or RRuntime()
    rscript = Check("Rscript", "Scientific analyses", runtime.available,
        runtime.version() or "" if runtime.available else "",
        "" if runtime.available else "Rscript could not be found.")
    packages = _r_packages(runtime.executable) if runtime.available else tuple(
        Check(name, R_REQUIRED_BY.get(name, R_OPTIONAL_BY.get(name, "")), False,
              detail="Rscript unavailable") for name in (*R_REQUIRED_BY, *R_OPTIONAL_BY))
    try:
        folder = r_scripts_dir()
        resources = tuple(Check(name, "Scientific execution", (folder / name).is_file(),
            detail="" if (folder / name).is_file() else "Required script missing")
            for name in CRITICAL_SCRIPTS)
    except FileNotFoundError as error:
        resources = tuple(Check(name, "Scientific execution", False, detail=str(error))
            for name in CRITICAL_SCRIPTS)
    writable = [_writable(Path(tempfile.gettempdir()), "Temporary directory")]
    if project_directory is not None:
        writable.append(_writable(Path(project_directory), "Project directory"))
    return Preflight(python, rscript, packages, resources, tuple(writable))
