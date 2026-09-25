from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .organism import get_organism
from .project import Project


SUPPORTED_ORGDB = {"9606": "org.Hs.eg.db", "10090": "org.Mm.eg.db"}
EVIDENCE_FILTERS = {"all", "exclude_iea", "experimental"}


@dataclass(frozen=True)
class GOReadiness:
    ready: bool
    reason: str


@dataclass(frozen=True)
class GOOutputs:
    annotations: pd.DataFrame
    unannotated: pd.DataFrame
    summary: pd.DataFrame
    metadata: dict[str, Any]
    tables: tuple[Path, ...]
    graphs: tuple[Path, ...]


def go_readiness(project: Project | None) -> GOReadiness:
    if project is None:
        return GOReadiness(False, "Open a project.")
    organism = get_organism(project)
    if organism is None:
        return GOReadiness(False, "Configure the project organism.")
    if organism.tax_id not in SUPPORTED_ORGDB:
        return GOReadiness(False, "Local GO is not configured for this organism yet.")
    if not (project.root / "mapping" / "tables" / "protein_catalog.csv").is_file():
        return GOReadiness(False, "Run mapping to obtain IDs usable by GO.")
    return GOReadiness(True, "Ready for local GO analysis.")


def _catalog(project: Project) -> pd.DataFrame:
    return pd.read_csv(project.root / "mapping" / "tables" / "protein_catalog.csv")


def available_sets(project: Project) -> dict[str, str]:
    sets = {"all_experiment": "All experiment proteins", "mapped": "Mapped proteins"}
    classification_path = project.root / "analyses" / "presence_absence" / "tables" / "classification.csv"
    if classification_path.is_file():
        frame = pd.read_csv(classification_path)
        for value in sorted(frame.get("classification", pd.Series(dtype=str)).dropna().unique()):
            sets[f"class:{value}"] = str(value)
    sets["manual"] = "Manual selection"
    return sets


def select_set(project: Project, selection: str, manual_rows: list[int] | None = None) -> pd.DataFrame:
    catalog = _catalog(project)
    if selection == "mapped":
        return catalog[catalog["mapping_status"] != "unmapped"].copy()
    if selection == "all_experiment":
        return catalog.copy()
    if selection == "manual":
        chosen = set(manual_rows or [])
        return catalog[catalog["source_row"].isin(chosen)].copy()
    if selection.startswith("class:"):
        classification = pd.read_csv(project.root / "analyses" / "presence_absence" / "tables" / "classification.csv")
        rows = classification.loc[classification["classification"] == selection[6:], "source_row"]
        return catalog[catalog["source_row"].isin(rows)].copy()
    raise ValueError("Unknown GO set.")


def entity_keys(frame: pd.DataFrame) -> set[str]:
    keys: set[str] = set()
    for column in ("ncbi_gene_id", "gene_symbol", "uniprot_accession"):
        if column in frame:
            keys.update(str(value) for value in frame[column].dropna() if str(value).strip())
    return keys


def prepare_go_arguments(project: Project, *, target_selection: str, background_selection: str,
    manual_rows: list[int] | None, ontologies: list[str], evidence_filter: str,
    fdr_cutoff: float, p_cutoff: float, min_count: int, top_n: int,
    simplify: bool, simplify_cutoff: float, run_id: str,
    allow_target_outside_background: bool = False) -> list[str]:
    state = go_readiness(project)
    if not state.ready:
        raise ValueError(state.reason)
    if not ontologies or any(item not in {"BP", "MF", "CC"} for item in ontologies):
        raise ValueError("Select BP, MF, and/or CC.")
    if evidence_filter not in EVIDENCE_FILTERS:
        raise ValueError("Invalid evidence filter.")
    target = select_set(project, target_selection, manual_rows)
    background = select_set(project, background_selection, manual_rows)
    target_keys, background_keys = entity_keys(target), entity_keys(background)
    outside = target_keys - background_keys
    if outside and not allow_target_outside_background:
        raise ValueError(f"{len(outside)} set identifier(s) are outside the selected background.")
    if outside:
        target = target[target.apply(lambda row: bool(entity_keys(pd.DataFrame([row])) & background_keys), axis=1)]
    raw = project.root / "analyses" / "GO" / "raw" / run_id
    if raw.exists():
        raise ValueError("GO run_id already exists.")
    raw.mkdir(parents=True)
    target_file, background_file = raw / "target.csv", raw / "background.csv"
    target.to_csv(target_file, index=False); background.to_csv(background_file, index=False)
    organism = get_organism(project); assert organism is not None
    return ["--input", str(project.root / project.config["input"]["processed_file"]),
        "--catalog", str(project.root / "mapping" / "tables" / "protein_catalog.csv"),
        "--output", str(project.root / "analyses" / "GO"), "--project", str(project.root),
        "--target-file", str(target_file), "--background-file", str(background_file),
        "--organism", organism.name, "--tax-id", organism.tax_id,
        "--ontologies", json.dumps(ontologies), "--evidence-filter", evidence_filter,
        "--fdr-cutoff", str(fdr_cutoff), "--p-cutoff", str(p_cutoff), "--min-count", str(min_count),
        "--top-n", str(top_n), "--simplify", str(simplify).lower(),
        "--simplify-cutoff", str(simplify_cutoff), "--run-id", run_id,
        "--target-name", target_selection, "--background-name", background_selection,
        "--outside-count", str(len(outside))]


def read_go_outputs(project: Project) -> GOOutputs:
    root = project.root / "analyses" / "GO"
    try:
        return GOOutputs(pd.read_csv(root / "annotation" / "go_annotations.csv"),
            pd.read_csv(root / "annotation" / "unannotated_proteins.csv"),
            pd.read_csv(root / "summary.csv"),
            json.loads((root / "latest_metadata.json").read_text(encoding="utf-8")),
            tuple(sorted(root.glob("**/*.csv"))), tuple(sorted((root / "graphs").glob("*.png"))))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Invalid GO results: {error}") from error


def list_go_runs(project: Project) -> list[str]:
    root = project.root / "analyses" / "GO" / "runs"
    return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []


def export_go(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, destination); return destination
