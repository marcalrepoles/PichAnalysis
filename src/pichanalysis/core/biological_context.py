"""Offline biological context for one frozen analysis run."""
from __future__ import annotations

import csv
import json
import re
import sqlite3
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .cross_module_integration import default_registry, open_cross_module_index
from .database_manager import DatabaseManager
from .entity_identity import CrossModuleMatch, EntityContext, MatchStatus
from . import differential_analysis, kegg_analysis, reactome_analysis, string_analysis

GREEN, RED, GRAY = "#258447", "#c43c3c", "#858b90"


class BiologicalContextError(RuntimeError):
    pass


@dataclass
class ContextItem:
    module: str
    run_id: str
    item_id: str
    title: str
    snapshot: str
    selected_entity: str
    members: pd.DataFrame
    edges: pd.DataFrame
    official_image: Path | None = None
    network_type: str = ""


def _text(value):
    value = str(value).strip()
    return "" if value.casefold() in {"", "nan", "none", "null", "na"} else value


def resolve_entity(project, module, run_id, identifier, source_row=None):
    if not _text(identifier):
        raise BiologicalContextError("Select a valid persisted entity.")
    index = open_cross_module_index(project)
    status = index.run_status(module, run_id)
    if not status or status[0] != "Ready":
        raise BiologicalContextError(f"Source run {module}/{run_id} is incomplete or unavailable.")
    records = [row for row in index.search(identifier)
               if row["module"] == module and row["run_id"] == run_id]
    if not records and module in {"kegg", "reactome", "string"}:
        key = {"kegg": "kegg_gene_id", "reactome": "reactome_entity_key",
               "string": "string_protein_id"}[module]
        records = [row for row in index.search(f"{key}:{identifier}")
                   if row["module"] == module and row["run_id"] == run_id]
    if source_row is not None:
        records = [row for row in records if row["source_row"] == int(source_row)]
    preferred = {"differential": "all_results", "mapping": "catalog",
                 "kegg": "mapping", "reactome": "mapping", "string": "mapping"}.get(module)
    if any(row["record_type"] == preferred for row in records):
        records = [row for row in records if row["record_type"] == preferred]
    if len(records) != 1:
        raise BiologicalContextError("No unique frozen source record was found. Select a specific source row.")
    return index, index.context(records[0]["id"])


def _details(index, project, match):
    record_id = match.details.get("record_id")
    if not record_id:
        return {}
    with sqlite3.connect(index.path) as db:
        row = db.execute("SELECT record_type,record_index FROM records WHERE id=? AND module=? AND run_id=?",
            (record_id, match.target_module, match.target_run_id)).fetchone()
    return default_registry().get(match.target_module).load_detail(
        project, match.target_run_id, row[0], row[1]) if row else {}


def _candidate_matches(index, context, include_other_lineages=False):
    matches = index.matches(context, include_other_lineages=include_other_lineages)
    if context.source_module in {"kegg", "reactome", "string"}:
        with sqlite3.connect(index.path) as db:
            db.row_factory = sqlite3.Row
            rows = db.execute("SELECT * FROM records WHERE module=? AND run_id=? AND record_type='mapping'",
                (context.source_module, context.source_run_id)).fetchall()
        own_ids = set()
        field = {"kegg": "kegg_gene_id", "reactome": "reactome_entity_key",
                 "string": "string_protein_id"}[context.source_module]
        for kind, value in context.module_specific_ids:
            if kind == field:
                own_ids.update(record["id"] for record in index.search(f"{kind}:{value}")
                    if record["module"] == context.source_module and
                    record["run_id"] == context.source_run_id)
        for row in rows:
            if (row["id"] in own_ids or
                    (context.source_row is not None and row["source_row"] == context.source_row)
                    or (context.feature_id and row["feature_id"] == context.feature_id)):
                matches.append(CrossModuleMatch(context.source_module, context.source_run_id,
                    row["target"] or "", MatchStatus.EXACT, "Source run record",
                    details={"record_id": str(row["id"])}))
    return [match for match in matches
            if match.target_module in {"kegg", "reactome", "string"}
            and match.status in {MatchStatus.EXACT, MatchStatus.MAPPED}]


def _official(module, snapshot, item_id, database_root=None):
    if not snapshot or not re.fullmatch(r"[A-Za-z0-9._-]+", snapshot):
        return None
    if not re.fullmatch(r"hsa[0-9]+" if module == "kegg" else r"R-HSA-[0-9]+", item_id):
        return None
    manager = DatabaseManager(database_root)
    if module == "kegg":
        path = manager.snapshots_root / snapshot / "images" / f"{item_id}.png"
    else:
        root = manager.reactome.snapshots / snapshot
        index = root / "diagram_index.csv"
        if not index.is_file():
            return None
        with index.open(newline="", encoding="utf-8") as stream:
            relative = next((row.get("relative_path", "") for row in csv.DictReader(stream)
                if row.get("Reactome_ID", "").upper() == item_id.upper()), "")
        base = (root / "diagrams").resolve()
        path = (base / relative).resolve() if relative else base
        if base not in path.parents:
            return None
    return path if path.is_file() else None


def find_contexts(project, index, context: EntityContext, database_root=None,
                  include_other_lineages=False):
    """Read membership and edges only from explicitly matched persisted runs."""
    grouped = {}
    for match in _candidate_matches(index, context, include_other_lineages):
        detail = _details(index, project, match)
        column = {"kegg": "kegg_gene_id", "reactome": "reactome_entity_key",
                  "string": "string_protein_id"}[match.target_module]
        selected = _text(detail.get(column))
        if selected:
            grouped.setdefault((match.target_module, match.target_run_id), set()).add(selected)
    items = []
    for (module, run_id), identities in sorted(grouped.items()):
        if len(identities) != 1:
            continue  # Never choose the first ambiguous biological identity.
        selected = next(iter(identities))
        status = index.run_status(module, run_id)
        if not status or status[0] != "Ready":
            raise BiologicalContextError(f"Target run {module}/{run_id} is incomplete.")
        if module == "string":
            output = string_analysis.read_string_outputs(project, run_id)
            nodes = output["tables"]["expanded_nodes"].fillna("").astype(str)
            edges = output["tables"]["expanded_edges"].fillna("").astype(str)
            if selected not in set(nodes.string_protein_id):
                continue
            direct = edges[(edges.protein_a == selected) | (edges.protein_b == selected)].copy()
            ids = {selected} | set(direct.protein_a) | set(direct.protein_b)
            members = nodes[nodes.string_protein_id.isin(ids)].copy()
            items.append(ContextItem(module, run_id, selected, f"STRING neighborhood: {selected}",
                _text(output["metadata"].get("snapshot_id")), selected, members, direct,
                network_type=_text(output["metadata"].get("network_type"))))
            continue
        if module == "kegg":
            output = kegg_analysis.read_kegg_outputs(project, run_id)
            path = output.root / "mapping/pathway_membership.csv"
            if not path.is_file():
                raise BiologicalContextError(f"Incomplete KEGG run {run_id}: pathway membership is missing.")
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            pathway, member, title = "pathway_id", "kegg_gene_id", "name"
            snapshot = _text(output.metadata.get("snapshot_id"))
        else:
            output = reactome_analysis.read_reactome_outputs(project, run_id)
            frame = output.membership.fillna("").astype(str)
            pathway, member, title = "Reactome_ID", "reactome_entity_key", "Pathway"
            snapshot = _text(output.metadata.get("snapshot_id"))
        if not {pathway, member}.issubset(frame.columns):
            raise BiologicalContextError(f"Incomplete {module} run {run_id}: membership columns are missing.")
        selected_rows = frame[frame[member] == selected]
        for item_id in sorted(selected_rows[pathway].dropna().astype(str).unique()):
            members = frame[frame[pathway] == item_id].copy()
            name = _text(members.iloc[0].get(title, item_id))
            if (module == "reactome" and hasattr(output, "hierarchy") and
                    {"child_pathway_id", "parent_pathway_name"}.issubset(output.hierarchy.columns)):
                parents = output.hierarchy[
                    output.hierarchy.child_pathway_id.astype(str) == item_id]
                parent_names = sorted(set(parents.parent_pathway_name.dropna().astype(str)) - {""})
                if parent_names:
                    name += " | Parent: " + ", ".join(parent_names[:3])
            items.append(ContextItem(module, run_id, item_id, name or item_id, snapshot,
                selected, members, pd.DataFrame(), _official(module, snapshot, item_id, database_root)))
    return items


def differential_choices(index, context):
    """Only matching frozen input hashes are automatically compatible."""
    source_hash = context.input_lineage.frozen_input_hash if context.input_lineage else None
    with sqlite3.connect(index.path) as db:
        rows = db.execute("SELECT run_id,status,input_hash FROM runs WHERE module='differential' ORDER BY run_id").fetchall()
    return [run_id for run_id, status, fingerprint in rows
            if status == "Ready" and ((context.source_module == "differential" and
                context.source_run_id == run_id) or (source_hash and fingerprint == source_hash))]



def other_differential_runs(index, compatible):
    """Offer unverified lineages only as explicit, caveated user choices."""
    with sqlite3.connect(index.path) as db:
        rows = db.execute("SELECT run_id,status FROM runs WHERE module='differential' ORDER BY run_id").fetchall()
    return [run_id for run_id, status in rows if status == "Ready" and run_id not in compatible]

def frozen_direction(project, run_id, context, member, differential_results=None):
    if not run_id:
        return "none"
    frame = (differential_results if differential_results is not None else
        differential_analysis.load_run(project, run_id)["tables"]["all_results"])
    frame = frame.fillna("").astype(str)
    selected = pd.DataFrame()
    member_uniprot = _text(member.get("UniProt", member.get("uniprot_accession", "")))
    member_gene = _text(member.get("Gene_symbol", member.get("gene_symbol", "")))
    same_selected_identity = (
        (member_uniprot and len(context.uniprot_accessions) == 1 and
         member_uniprot == context.uniprot_accessions[0]) or
        (not member_uniprot and member_gene and len(context.gene_symbols) == 1 and
         member_gene == context.gene_symbols[0]))
    if (context.source_module == "differential" and context.source_run_id == run_id
            and same_selected_identity):
        if context.source_row is not None and "source_row" in frame:
            selected = frame[frame.source_row == str(context.source_row)]
        elif context.feature_id and "feature_id" in frame:
            selected = frame[frame.feature_id == context.feature_id]
    if selected.empty:
        for column, value in (("uniprot_accession", member_uniprot), ("gene_symbol", member_gene)):
            if value and column in frame:
                candidate = frame[frame[column] == value]
                if len(candidate) == 1:
                    selected = candidate
                    break
    if len(selected) != 1:
        return "none"
    row = selected.iloc[0]
    if _text(row.get("combined_significant")).lower() != "true":
        return "none"
    return {"higher_in_condition_A": "up", "higher_in_condition_B": "down"}.get(
        _text(row.get("effect_direction")), "none")


def save_context_figure(project, context, item, differential_run, image):
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "-" + uuid.uuid4().hex[:8]
    root = Path(project.root) / "analyses/biological_context/runs" / run_id
    for folder in ("tables", "figures", "provenance"):
        (root / folder).mkdir(parents=True, exist_ok=False)
    item.members.to_csv(root / "tables/members.csv", index=False)
    item.edges.to_csv(root / "tables/edges.csv", index=False)
    figure = root / "figures/context.png"
    if not image.save(str(figure), "PNG"):
        raise BiologicalContextError("Could not save the generated context figure.")
    metadata = {"source_entity": context.original_identifier or context.feature_id or "",
        "source_module": context.source_module, "source_run": context.source_run_id,
        "source_row": context.source_row, "differential_run": differential_run,
        "item_module": item.module, "item_run": item.run_id, "item_id": item.item_id,
        "snapshot": item.snapshot, "generation_timestamp": datetime.now(timezone.utc).isoformat(),
        "legend": {"up": GREEN, "down": RED, "no_direction": GRAY},
        "interpretation": "Status of detected members in this sample/run; not pathway activity.",
        "figure": "figures/context.png"}
    (root / "provenance/context.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return root


def list_saved_contexts(project):
    root = Path(project.root) / "analyses/biological_context/runs"
    return sorted((path for path in root.iterdir() if (path / "provenance/context.json").is_file()),
        reverse=True) if root.is_dir() else []


def load_saved_context(root):
    root = Path(root)
    metadata = json.loads((root / "provenance/context.json").read_text(encoding="utf-8"))
    figure = root / metadata["figure"]
    if not figure.is_file() or not (root / "tables/members.csv").is_file():
        raise BiologicalContextError("Historical biological context is incomplete.")
    return metadata, figure