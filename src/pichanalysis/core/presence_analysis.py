from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .column_mapping import ColumnRole, experimental_design
from .project import Project


@dataclass(frozen=True)
class PresenceReadiness:
    ready: bool
    reason: str


@dataclass(frozen=True)
class PresenceOutputs:
    presence_matrix: pd.DataFrame
    condition_detection: pd.DataFrame
    classification: pd.DataFrame
    condition_summary: pd.DataFrame
    metadata: dict[str, Any]
    graphs: tuple[Path, ...]


def quantitative_columns(project: Project) -> list[dict[str, str]]:
    return [
        {"column": name, "condition": str(item.get("condition", "")),
         "replicate": str(item.get("replicate", "")),
         "quantification_type": str(item.get("quantification_type", "unknown"))}
        for name, item in project.config.get("columns", {}).items()
        if item.get("role") == ColumnRole.QUANTIFICATION.value
    ]


def quantification_types(project: Project) -> list[str]:
    return sorted({item["quantification_type"] for item in quantitative_columns(project)})


def presence_readiness(project: Project | None, quantification_type: str | None = None) -> PresenceReadiness:
    if project is None:
        return PresenceReadiness(False, "Abra um projeto.")
    design = experimental_design(project)
    if not design["primary_identifier_column"]:
        return PresenceReadiness(False, "Configure um identificador principal.")
    columns = quantitative_columns(project)
    if not columns:
        return PresenceReadiness(False, "Configure ao menos uma coluna Quantification.")
    if any(not item["condition"] for item in columns):
        return PresenceReadiness(False, "Todas as colunas quantitativas precisam de uma condição.")
    selected = quantification_type or (quantification_types(project)[0] if len(quantification_types(project)) == 1 else None)
    if not selected:
        return PresenceReadiness(False, "Escolha o tipo de quantificação.")
    chosen = [item for item in columns if item["quantification_type"] == selected]
    if not chosen:
        return PresenceReadiness(False, "O tipo de quantificação selecionado não possui colunas.")
    if not project.config.get("input", {}).get("processed_file"):
        return PresenceReadiness(False, "Importe uma tabela antes da análise.")
    return PresenceReadiness(True, "Pronto para executar presença/ausência.")


def replicate_counts(project: Project, quantification_type: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in quantitative_columns(project):
        if item["quantification_type"] == quantification_type:
            counts[item["condition"]] = counts.get(item["condition"], 0) + 1
    return counts


def suggested_minimum(project: Project, quantification_type: str) -> int:
    counts = replicate_counts(project, quantification_type)
    return max(1, math.ceil(min(counts.values()) / 2)) if counts else 1


def build_presence_arguments(project: Project, *, quantification_type: str,
    selected_conditions: list[str], zero_is_missing: bool, threshold: float,
    rule_mode: str, rule_value: float, predominant: bool, run_id: str) -> list[str]:
    state = presence_readiness(project, quantification_type)
    if not state.ready:
        raise ValueError(state.reason)
    available = replicate_counts(project, quantification_type)
    if not selected_conditions or any(item not in available for item in selected_conditions):
        raise ValueError("Selecione ao menos uma condição válida.")
    if rule_mode not in {"count", "fraction"}:
        raise ValueError("Modo de reprodutibilidade inválido.")
    if rule_value <= 0 or (rule_mode == "fraction" and rule_value > 1):
        raise ValueError("Critério de reprodutibilidade inválido.")
    design = experimental_design(project)
    column_map = [item for item in quantitative_columns(project)
                  if item["quantification_type"] == quantification_type and item["condition"] in selected_conditions]
    return ["--input", str(project.root / project.config["input"]["processed_file"]),
        "--output", str(project.root / "analyses" / "presence_absence"),
        "--project", str(project.root), "--identifier-column", str(design["primary_identifier_column"]),
        "--identifier-type", str(design["primary_identifier_type"]),
        "--quantification-type", quantification_type,
        "--column-map", json.dumps(column_map, ensure_ascii=False, separators=(",", ":")),
        "--conditions", json.dumps(selected_conditions, ensure_ascii=False, separators=(",", ":")),
        "--zero-is-missing", str(zero_is_missing).lower(), "--threshold", str(float(threshold)),
        "--rule-mode", rule_mode, "--rule-value", str(float(rule_value)),
        "--predominant", str(predominant).lower(), "--run-id", run_id]


def read_presence_outputs(project: Project) -> PresenceOutputs:
    root = project.root / "analyses" / "presence_absence"
    try:
        metadata = json.loads((root / "latest_metadata.json").read_text(encoding="utf-8"))
        tables = root / "tables"
        return PresenceOutputs(pd.read_csv(tables / "presence_matrix.csv"),
            pd.read_csv(tables / "condition_detection.csv"),
            pd.read_csv(tables / "classification.csv"),
            pd.read_csv(tables / "condition_summary.csv"), metadata,
            tuple(sorted((root / "graphs").glob("*.png"))))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Resultados de presença/ausência inválidos: {error}") from error


def list_presence_runs(project: Project) -> list[str]:
    root = project.root / "analyses" / "presence_absence" / "runs"
    return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []


def export_presence(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination
