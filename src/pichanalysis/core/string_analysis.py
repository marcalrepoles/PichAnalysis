"""Offline STRING run preparation and persisted-result loading.

Scientific graph measures are computed by r_scripts/08_string_analysis.R.
"""
from __future__ import annotations

import csv
import json
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .r_runtime import RRuntime


PRESETS = {"Exploratory": 150, "Functional": 400, "Robust": 700}
TABLES = {
    "mapping": "mapping/string_mapping.csv", "ambiguous": "mapping/ambiguous.csv",
    "unmapped": "mapping/unmapped.csv", "seeds": "networks/seeds.csv",
    "internal_nodes": "networks/internal_nodes.csv", "internal_edges": "networks/internal_edges.csv",
    "degree1_nodes": "networks/degree1_nodes.csv", "degree1_edges": "networks/degree1_edges.csv",
    "degree2_nodes": "networks/degree2_nodes.csv", "degree2_edges": "networks/degree2_edges.csv",
    "expanded_nodes": "networks/expanded_nodes.csv", "expanded_edges": "networks/expanded_edges.csv",
    "common_direct_neighbors": "networks/common_direct_neighbors.csv",
    "seed_support": "networks/seed_support.csv", "seed_seed_edges": "networks/seed_seed_edges.csv",
    "edge_evidence": "networks/edge_evidence.csv", "node_metrics": "metrics/node_metrics.csv",
    "internal_node_metrics": "metrics/internal_node_metrics.csv", "internal_components": "metrics/internal_components.csv",
    "seed_metrics": "metrics/seed_metrics.csv", "hubs": "metrics/hubs.csv",
    "components": "metrics/components.csv", "evidence_summary": "metrics/evidence_summary.csv",
    "network_summary": "metrics/network_summary.csv", "summary": "summary.csv",
}
PLOTS = ("top_degree", "top_betweenness", "component_sizes", "degree_distribution",
         "seed_direct_neighbors", "seed_support_distribution")


class StringAnalysisError(RuntimeError): pass
class StringDatabaseUnavailableError(StringAnalysisError): pass
class UnsupportedOrganismError(StringAnalysisError): pass
class NoStringSeedsError(StringAnalysisError): pass
class StringPresenceRequiredError(StringAnalysisError): pass
class StringTargetOutsideBackgroundError(StringAnalysisError):
    def __init__(self, outside, target_size, background_size):
        self.details = {"entities_outside_background": sorted(outside), "target_size": target_size,
                        "background_size": background_size}
        super().__init__(f"Target outside background: {len(outside)} protein(s) require authorization.")
class StringExpansionLimitError(StringAnalysisError):
    def __init__(self, observed, limit, threshold, hop, network_type):
        self.details = dict(observed=observed, limit=limit, threshold=threshold, hop=hop,
                            network_type=network_type)
        super().__init__(f"External expansion exceeds safety limit: {observed} > {limit} nodes at hop {hop}.")
class StringIgraphUnavailableError(StringAnalysisError): pass
class StringRExecutionError(StringAnalysisError): pass
class MissingStringOutputError(StringAnalysisError): pass


@dataclass(frozen=True)
class StringParameters:
    target_selection: str = "All mapped proteins"
    background_selection: str = "All mapped proteins"
    manual_rows: tuple[int, ...] = ()
    background_manual_rows: tuple[int, ...] = ()
    network_type: str = "functional"
    threshold_preset: str = "Functional"
    combined_score_threshold: int | None = None
    max_hop: int = 1
    degree1_selection_mode: str = "union"
    max_external_nodes: int = 10000
    minimum_hub_degree: int = 3
    top_n: int = 20
    allow_target_outside_background: bool = False

    @property
    def threshold(self):
        if self.threshold_preset == "Custom":
            if self.combined_score_threshold is None: raise StringAnalysisError("Custom threshold requires an integer score.")
            return self.combined_score_threshold
        if self.threshold_preset not in PRESETS: raise StringAnalysisError("Unknown PichAnalysis threshold preset.")
        return PRESETS[self.threshold_preset]


def _root(project): return Path(project.root)
def _catalog(project):
    path = _root(project) / "mapping/tables/protein_catalog.csv"
    if not path.is_file(): raise StringAnalysisError("Experimental protein catalog is unavailable.")
    return pd.read_csv(path, dtype=str).fillna("")
def _tokens(value):
    import re
    return sorted({x for x in re.split(r"[;,|\s]+", str(value).strip()) if x and x.upper() != "NA"})
def _rows_for_selection(project, catalog, selection, manual_rows):
    if selection in {"All mapped proteins", "all_mapped"}: return catalog
    if selection in {"Manual selection", "manual"}: return catalog[catalog.source_row.isin({str(x) for x in manual_rows})]
    path = _root(project) / "analyses/presence_absence/tables/classification.csv"
    if not path.is_file(): raise StringPresenceRequiredError("Presence/Absence outputs are required for this target.")
    classes = pd.read_csv(path, dtype=str).fillna("")
    label = selection[6:] if selection.startswith("class:") else selection
    if selection == "Reproducibly detected":
        selected = classes.loc[classes.classification != "Not reproducibly detected", "source_row"]
    else: selected = classes.loc[classes.classification == label, "source_row"]
    return catalog[catalog.source_row.isin(set(selected))]


def map_experiment(project, database):
    catalog = _catalog(project); records = []
    for source_row, group in catalog.groupby("source_row", sort=False):
        accessions = sorted({accession for value in group.get("uniprot_accession", pd.Series(dtype=str)) for accession in _tokens(value)})
        candidates = {}
        for accession in accessions:
            lookup = database.lookup_uniprot(accession)
            for string_id in lookup["candidate_string_ids"]:
                candidates.setdefault(string_id, set()).update(lookup["candidate_sources"])
        status = "mapped_unique" if len(candidates) == 1 else "ambiguous" if candidates else "unmapped"
        first = group.iloc[0]
        for string_id in (sorted(candidates) or [""]):
            protein = database.get_protein(string_id) if string_id else None
            records.append({"source_row": source_row, "original_id": first.get("original_id", ""),
                            "uniprot_accession": ";".join(accessions), "gene_symbol": first.get("gene_symbol", ""),
                            "canonical_protein_key": "UP:" + accessions[0] if len(accessions) == 1 else "",
                            "string_protein_id": string_id, "preferred_name": (protein or {}).get("preferred_name", ""),
                            "mapping_status": status, "mapping_source": ";".join(sorted(candidates.get(string_id, ()))),
                            "candidate_string_ids": ";".join(sorted(candidates)), "source_rows": source_row})
    columns = ("source_row original_id uniprot_accession gene_symbol canonical_protein_key string_protein_id "
               "preferred_name mapping_status mapping_source candidate_string_ids source_rows").split()
    frame = pd.DataFrame(records, columns=columns)
    for string_id, rows in frame.loc[frame.mapping_status == "mapped_unique"].groupby("string_protein_id"):
        frame.loc[rows.index, "source_rows"] = ";".join(sorted(set(rows.source_row), key=str))
    return frame


def _edge_columns(conn, table): return [row[1] for row in conn.execute(f"PRAGMA table_info({table})")]
def _incident(conn, table, ids, threshold):
    """Use the endpoint indexes; never scan or load the whole STRING graph."""
    if not ids: return []
    rows = {}
    for node in ids:
        for row in conn.execute(f"SELECT * FROM {table} WHERE protein_a=? AND combined_score>=? UNION ALL "
                                f"SELECT * FROM {table} WHERE protein_b=? AND combined_score>=?",
                                (node, threshold, node, threshold)):
            rows[(row["protein_a"], row["protein_b"])] = dict(row)
    return list(rows.values())
def _neighbor(edge, node): return edge["protein_b"] if edge["protein_a"] == node else edge["protein_a"]
def _write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def prepare_string_run(project, manager, run_id, parameters=None):
    p = parameters or StringParameters(); root = _root(project)
    if p.network_type not in {"functional", "physical"}: raise StringAnalysisError("Network type must be functional or physical.")
    threshold = p.threshold
    if not isinstance(threshold, int) or not 0 <= threshold <= 1000: raise StringAnalysisError("Combined score threshold must be an integer from 0 to 1000.")
    if p.max_hop not in (0, 1, 2): raise StringAnalysisError("Maximum hop must be 0, 1 or 2.")
    if p.degree1_selection_mode not in {"union", "strict_common"}: raise StringAnalysisError("Unknown Degree 1 selection mode.")
    if p.max_external_nodes < 0 or p.minimum_hub_degree < 0 or p.top_n < 1: raise StringAnalysisError("Invalid operational or plot limit.")
    snapshot = manager.string.active_snapshot()
    if snapshot is None: raise StringDatabaseUnavailableError("STRING database unavailable.")
    manifest = manager.string.manifest(snapshot)
    if str(manifest.get("tax_id")) != "9606": raise UnsupportedOrganismError("Only Homo sapiens (9606) is supported.")
    if str(getattr(project, "config", {}).get("organism_tax_id") or "9606") != "9606":
        raise UnsupportedOrganismError("Only Homo sapiens (9606) is supported.")
    catalog = _catalog(project); mapping = map_experiment(project, manager.string)
    target_rows = set(_rows_for_selection(project, catalog, p.target_selection, p.manual_rows).source_row)
    background_rows = set(_rows_for_selection(project, catalog, p.background_selection, p.background_manual_rows).source_row)
    target_ids = set(mapping.loc[mapping.source_row.isin(target_rows) & (mapping.mapping_status == "mapped_unique"), "string_protein_id"])
    background_ids = set(mapping.loc[mapping.source_row.isin(background_rows) & (mapping.mapping_status == "mapped_unique"), "string_protein_id"])
    outside = target_ids - background_ids
    if outside and not p.allow_target_outside_background:
        raise StringTargetOutsideBackgroundError(outside, len(target_ids), len(background_ids))
    seeds = sorted(target_ids & background_ids)
    if not seeds: raise NoStringSeedsError("No uniquely mapped STRING seeds are available.")
    base = root / "analyses/STRING"; run = base / "runs" / run_id
    provenance = root / "scripts/runs" / f"{run_id}_string"
    if run.exists() or provenance.exists(): raise StringAnalysisError("This STRING run ID already exists.")
    table = p.network_type + "_edges"; db_path = snapshot / "string_network.sqlite"
    with sqlite3.connect(f"file:{db_path.as_posix()}?mode=ro", uri=True) as conn:
        conn.row_factory = sqlite3.Row; edge_columns = _edge_columns(conn, table)
        first_edges = _incident(conn, table, seeds, threshold)
        direct = {}
        seed_set = set(seeds)
        for edge in first_edges:
            for seed in (edge["protein_a"], edge["protein_b"]):
                if seed in seed_set:
                    neighbor = _neighbor(edge, seed)
                    if neighbor not in seed_set: direct.setdefault(neighbor, {})[seed] = edge["combined_score"]
        common = {node: support for node, support in direct.items() if len(support) == len(seeds)}
        degree1 = sorted(common if p.degree1_selection_mode == "strict_common" else direct) if p.max_hop >= 1 else []
        if len(degree1) > p.max_external_nodes: raise StringExpansionLimitError(len(degree1), p.max_external_nodes, threshold, 1, p.network_type)
        second_edges = _incident(conn, table, degree1, threshold) if p.max_hop >= 2 else []
        degree2 = sorted({endpoint for edge in second_edges for endpoint in (edge["protein_a"], edge["protein_b"])} - set(seeds) - set(degree1))
        if len(degree1) + len(degree2) > p.max_external_nodes:
            raise StringExpansionLimitError(len(degree1) + len(degree2), p.max_external_nodes, threshold, 2, p.network_type)
        included = set(seeds) | set(degree1) | set(degree2)
        all_edges = { (e["protein_a"], e["protein_b"]): e for e in first_edges + second_edges }
        # Add edges among included nodes, including neighbor-neighbor links, through indexed endpoint queries.
        if p.max_hop >= 1:
            all_edges.update({(e["protein_a"], e["protein_b"]): e for e in _incident(conn, table, included, threshold)})
        expanded_edges = sorted((e for e in all_edges.values() if e["protein_a"] in included and e["protein_b"] in included), key=lambda e:(e["protein_a"],e["protein_b"]))
        internal_edges = [e for e in expanded_edges if e["protein_a"] in seeds and e["protein_b"] in seeds]
        degree1_edges = [e for e in expanded_edges if (e["protein_a"] in degree1 or e["protein_b"] in degree1) and e["protein_a"] in set(seeds)|set(degree1) and e["protein_b"] in set(seeds)|set(degree1)]
        degree2_edges = [e for e in expanded_edges if e["protein_a"] in degree2 or e["protein_b"] in degree2]
        # Local metadata only. Temporary IN lists stay below SQLite's parameter limit.
        info = {}
        for node in sorted(included):
            row = conn.execute("SELECT * FROM proteins WHERE string_protein_id=?", (node,)).fetchone()
            info[node] = dict(row) if row else {}
    classification_path = root / "analyses/presence_absence/tables/classification.csv"
    classifications = {}
    if classification_path.is_file():
        classes = pd.read_csv(classification_path, dtype=str).fillna("")
        classifications = dict(zip(classes.source_row, classes.classification))
    seed_meta = {}
    for node in seeds:
        rows = mapping[(mapping.string_protein_id == node) & (mapping.mapping_status == "mapped_unique")]
        seed_meta[node] = {"uniprot_accession": ";".join(sorted({x for value in rows.uniprot_accession for x in _tokens(value)})),
                           "gene_symbol": ";".join(sorted(set(filter(None, rows.gene_symbol)))),
                           "canonical_protein_key": ";".join(sorted(set(filter(None, rows.canonical_protein_key)))),
                           "source_rows": ";".join(sorted({x for value in rows.source_rows for x in value.split(";")}))}
        seed_meta[node]["presence_absence_class"] = ";".join(sorted({classifications[x] for x in seed_meta[node]["source_rows"].split(";") if x in classifications and classifications[x]}))
    seed_support = []
    for node in degree1:
        for seed in sorted(direct[node]):
            seed_support.append({"string_protein_id": node, "hop_level": 1, "supporting_seed": seed_meta[seed]["uniprot_accession"],
                                 "supporting_seed_string_id": seed, "supporting_seed_gene": seed_meta[seed]["gene_symbol"]})
    # For hop 2, support means a selected seed can reach the node through a selected Degree 1 parent.
    parents = {}
    degree1_set, degree2_set = set(degree1), set(degree2)
    for edge in second_edges:
        for parent in (edge["protein_a"], edge["protein_b"]):
            if parent in degree1_set:
                child = _neighbor(edge, parent)
                if child in degree2_set: parents.setdefault(child, set()).update(direct[parent])
    for node in degree2:
        for seed in sorted(parents.get(node, ())):
            seed_support.append({"string_protein_id": node, "hop_level": 2, "supporting_seed": seed_meta[seed]["uniprot_accession"],
                                 "supporting_seed_string_id": seed, "supporting_seed_gene": seed_meta[seed]["gene_symbol"]})
    support_sets = {}
    for row in seed_support: support_sets.setdefault(row["string_protein_id"], set()).add(row["supporting_seed_string_id"])
    support_count = {node: len(support_sets.get(node, ())) for node in degree1+degree2}
    nodes = []
    for node in sorted(included):
        hop = 0 if node in seeds else 1 if node in degree1 else 2
        nodes.append({"string_protein_id": node, "preferred_name": info[node].get("preferred_name", ""),
                      "annotation": info[node].get("annotation", ""), "gene_symbol": seed_meta.get(node, {}).get("gene_symbol", "") or info[node].get("preferred_name", ""),
                      "uniprot_accession": seed_meta.get(node, {}).get("uniprot_accession", ""),
                      "presence_absence_class": seed_meta.get(node, {}).get("presence_absence_class", ""),
                      "canonical_protein_key": seed_meta.get(node, {}).get("canonical_protein_key", ""),
                      "source_rows": seed_meta.get(node, {}).get("source_rows", ""), "hop_level": hop,
                      "is_experimental_protein": node in seeds, "seed_support_count": support_count.get(node, 0)})
    node_columns = list(nodes[0]); support_columns = ["string_protein_id", "hop_level", "supporting_seed", "supporting_seed_string_id", "supporting_seed_gene"]
    common_rows = []
    for node, support in sorted(common.items()):
        values = list(support.values()); common_rows.append(dict(string_protein_id=node, preferred_name=info.get(node,{}).get("preferred_name",""), gene_symbol=info.get(node,{}).get("preferred_name",""), seed_count=len(support), required_seed_count=len(seeds), supporting_seeds=";".join(sorted(support)), minimum_edge_score_to_seed=min(values), maximum_edge_score_to_seed=max(values), mean_edge_score_to_seed=sum(values)/len(values)))
    common_columns = ["string_protein_id", "preferred_name", "gene_symbol", "seed_count", "required_seed_count", "supporting_seeds", "minimum_edge_score_to_seed", "maximum_edge_score_to_seed", "mean_edge_score_to_seed"]
    run.mkdir(parents=True); provenance.mkdir(parents=True)
    for name, rows, columns in (
        ("networks/seeds.csv", [x for x in nodes if x["hop_level"]==0], node_columns),
        ("networks/internal_nodes.csv", [x for x in nodes if x["hop_level"]==0], node_columns),
        ("networks/degree1_nodes.csv", [x for x in nodes if x["hop_level"]==1], node_columns),
        ("networks/degree2_nodes.csv", [x for x in nodes if x["hop_level"]==2], node_columns),
        ("networks/expanded_nodes.csv", nodes, node_columns),
        ("networks/internal_edges.csv", internal_edges, edge_columns),
        ("networks/degree1_edges.csv", degree1_edges, edge_columns),
        ("networks/degree2_edges.csv", degree2_edges, edge_columns),
        ("networks/expanded_edges.csv", expanded_edges, edge_columns),
        ("networks/seed_seed_edges.csv", internal_edges, edge_columns),
        ("networks/edge_evidence.csv", expanded_edges, edge_columns),
        ("networks/seed_support.csv", seed_support, support_columns),
        ("networks/common_direct_neighbors.csv", common_rows, common_columns),
    ): _write_csv(run/name, rows, columns)
    for name, frame in (("string_mapping", mapping), ("ambiguous", mapping[mapping.mapping_status=="ambiguous"]), ("unmapped", mapping[mapping.mapping_status=="unmapped"])):
        path = run/"mapping"/f"{name}.csv"; path.parent.mkdir(parents=True,exist_ok=True); frame.to_csv(path,index=False)
    metadata = dict(run_id=run_id, created_at=datetime.now(timezone.utc).isoformat(), snapshot_id=snapshot.name,
                    snapshot_path=str(snapshot), string_version=manifest["string_version"], sqlite_sha256=manifest["sqlite_sha256"],
                    raw_source_hashes={x["filename"]:x["sha256"] for x in manifest["files"]}, network_type=p.network_type,
                    combined_score_threshold=threshold, threshold_preset=p.threshold_preset, max_hop=p.max_hop,
                    degree1_selection_mode=p.degree1_selection_mode, target_definition=p.target_selection,
                    background_definition=p.background_selection, minimum_hub_degree=p.minimum_hub_degree,
                    max_external_nodes=p.max_external_nodes, top_n=p.top_n, seed_count=len(seeds),
                    degree1_count=len(degree1), degree2_count=len(degree2), metrics_graph_scope="internal" if p.max_hop==0 else f"expanded through hop {p.max_hop}",
                    methodology={"betweenness":"undirected unweighted topology", "closeness":"igraph component-local normalized closeness", "seed_support":"direct qualifying edge for hop 1; selected two-hop path for hop 2", "network_access":False})
    (run/"metadata.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    (provenance/"string_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    (provenance/"parameters.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    scripts = Path(__file__).resolve().parents[3]/"r_scripts"
    for source in (scripts/"08_string_analysis.R", scripts/"lib/string_analysis.R"):
        shutil.copy2(source,provenance/source.name)
    return run, provenance


def read_string_outputs(project, run_id):
    run = _root(project)/"analyses/STRING/runs"/run_id
    missing = [str(run/relative) for relative in TABLES.values() if not (run/relative).is_file()]
    missing += [str(run/name) for name in ("metadata.json","STRING_analysis.xlsx") if not (run/name).is_file()]

    if missing: raise MissingStringOutputError("Missing expected STRING output: "+", ".join(missing))
    return {"run_root":run, "metadata":json.loads((run/"metadata.json").read_text(encoding="utf-8")),
            "tables":{name:pd.read_csv(run/path) for name,path in TABLES.items()},
            "workbook":run/"STRING_analysis.xlsx", "plots":tuple(sorted((run/"plots").glob("*.*")))}


def list_string_runs(project):
    folder = _root(project)/"analyses/STRING/runs"
    return sorted(p.name for p in folder.iterdir() if p.is_dir()) if folder.is_dir() else []


def run_string_analysis(project, manager, runtime:RRuntime, *, run_id, parameters=None, script=None, timeout=600):
    run, provenance = prepare_string_run(project,manager,run_id,parameters)
    entry = script or Path(__file__).resolve().parents[3]/"r_scripts/08_string_analysis.R"
    try: result = runtime.run(entry,"--run",str(run),"--provenance",str(provenance),timeout=timeout)
    except RuntimeError as exc: raise StringRExecutionError("R execution failure: "+str(exc)) from exc
    if result.returncode:
        message=(result.stderr or result.stdout or "Rscript exited with an error.").strip()
        if "igraph package is required" in message: raise StringIgraphUnavailableError(message)
        raise StringRExecutionError("R execution failure: "+message)
    outputs = read_string_outputs(project,run_id)
    latest = _root(project) / "analyses/STRING"
    (latest / "mapping").mkdir(parents=True, exist_ok=True)
    for name in ("string_mapping.csv", "ambiguous.csv", "unmapped.csv"):
        shutil.copy2(run / "mapping" / name, latest / "mapping" / name)
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "STRING_analysis.xlsx", latest / "STRING_analysis.xlsx")
    return outputs


def string_node_neighbors(outputs, string_id):
    """Filter persisted run edges only; no SQLite or network access."""
    edges = outputs["tables"]["expanded_edges"]
    return edges[(edges.protein_a.astype(str)==str(string_id)) | (edges.protein_b.astype(str)==str(string_id))].copy()


def string_node_support(outputs, string_id):
    support = outputs["tables"]["seed_support"]
    return support[support.string_protein_id.astype(str)==str(string_id)].copy()


def string_edge_details(outputs, protein_a, protein_b):
    edges = outputs["tables"]["edge_evidence"]
    mask = ((edges.protein_a.astype(str)==str(protein_a)) & (edges.protein_b.astype(str)==str(protein_b))) | ((edges.protein_a.astype(str)==str(protein_b)) & (edges.protein_b.astype(str)==str(protein_a)))
    return edges[mask].copy()