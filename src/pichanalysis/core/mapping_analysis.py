from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from platformdirs import user_cache_path

from .column_mapping import experimental_design, validate_mapping
from .organism import get_organism
from .project import Project


@dataclass(frozen=True)
class AnalysisReadiness:
    ready: bool
    reason: str


@dataclass(frozen=True)
class MappingOutputs:
    catalog: pd.DataFrame
    mapping: pd.DataFrame
    unmapped: pd.DataFrame
    ambiguous: pd.DataFrame
    metadata: dict[str, Any]


def mapping_readiness(project: Project | None) -> AnalysisReadiness:
    if project is None:
        return AnalysisReadiness(False, "Abra um projeto.")
    if get_organism(project) is None:
        return AnalysisReadiness(False, "Configure o organismo do projeto.")
    columns = project.config.get("columns", {})
    validation = validate_mapping(columns if isinstance(columns, dict) else {})
    if not validation.valid:
        return AnalysisReadiness(False, validation.errors[0])
    design = experimental_design(project)
    if not design["primary_identifier_column"] or not design["primary_identifier_type"]:
        return AnalysisReadiness(False, "Configure um identificador principal.")
    if not project.config.get("input", {}).get("processed_file"):
        return AnalysisReadiness(False, "Importe uma tabela antes de mapear.")
    return AnalysisReadiness(True, "Pronto para mapear e anotar.")


def cache_path() -> Path:
    path = user_cache_path("PichAnalysis", "PichAnalysis") / "annotation_cache.sqlite"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def new_run_id() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def build_mapping_arguments(project: Project, cache: Path, refresh: bool, run_id: str) -> list[str]:
    readiness = mapping_readiness(project)
    if not readiness.ready:
        raise ValueError(readiness.reason)
    design = experimental_design(project)
    organism = get_organism(project)
    assert organism is not None
    input_path = project.root / project.config["input"]["processed_file"]
    return ["--input", str(input_path), "--output", str(project.root / "mapping"),
        "--id-column", str(design["primary_identifier_column"]),
        "--id-type", str(design["primary_identifier_type"]), "--tax-id", organism.tax_id,
        "--organism", organism.name, "--cache", str(cache),
        "--refresh", "true" if refresh else "false", "--run-id", run_id,
        "--project", str(project.root)]


def read_mapping_outputs(project: Project) -> MappingOutputs:
    tables = project.root / "mapping" / "tables"
    metadata_path = project.root / "mapping" / "latest_metadata.json"
    try:
        return MappingOutputs(pd.read_csv(tables / "protein_catalog.csv"),
            pd.read_csv(tables / "id_mapping.csv"), pd.read_csv(tables / "unmapped.csv"),
            pd.read_csv(tables / "ambiguous.csv"),
            json.loads(metadata_path.read_text(encoding="utf-8")))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Resultados de mapeamento ausentes ou inválidos: {error}") from error


def export_result(source: Path, destination: Path) -> Path:
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        raise FileNotFoundError(f"Resultado não encontrado: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination

