"""Conservative, offline identity and index behavior."""
import json
import sqlite3

import pandas as pd

from pichanalysis.core.cross_module_integration import (CrossModuleAdapter,
    CrossModuleAdapterRegistry, open_cross_module_index, rebuild_cross_module_index)
from pichanalysis.core.entity_identity import InputLineage, MatchStatus, tokens
from pichanalysis.core.project import create_project


def _registry(state):
    registry = CrossModuleAdapterRegistry()
    for module in ("alpha", "beta"):
        registry.register(CrossModuleAdapter(module,
            lambda project, name=module: list(state[name]),
            lambda project, run_id, name=module: {
                "run_root": project.root / "analyses" / name / "runs" / run_id,
                "metadata": {"run_id": run_id, "snapshot_id": state[name][run_id]["snapshot"]},
                "tables": {"records": state[name][run_id]["frame"]}},
            ("records",)))
    return registry


def _frame(*rows):
    return pd.DataFrame(rows, columns=["feature_id", "source_row", "original_identifier",
        "uniprot_accession", "gene_symbol", "ncbi_gene_id"])


def test_identity_types_and_isoforms_remain_distinct(tmp_path):
    project = create_project(tmp_path, "integration")
    state = {
        "alpha": {"source": {"snapshot": "X", "frame": _frame(
            ("ROW000001", "1", "P12345-2", "P12345-2", "GENE", "100"))}},
        "beta": {"target": {"snapshot": "Y", "frame": _frame(
            ("ROW000002", "2", "P12345", "P12345", "GENE", "101"),
            ("ROW000003", "3", "P12345-2", "P12345-2", "GENE", "102"))}},
    }
    registry = _registry(state)
    index = rebuild_cross_module_index(project, registry)
    source = next(row for row in index.search("P12345-2") if row["module"] == "alpha")
    context = index.context(source["id"])
    assert context.uniprot_accessions == ("P12345-2",)
    assert context.gene_symbols == ("GENE",)
    assert context.ncbi_gene_ids == ("100",)
    assert not index.matches(context)
    cross = index.matches(context, include_other_lineages=True)
    assert any(match.source_identifier == "P12345-2" for match in cross)
    assert not any(match.source_identifier == "P12345" for match in cross)
    assert all(not match.lineage_compatible for match in cross)
    assert {match.database_snapshot for match in cross} == {"Y"}
    assert tokens("P11111;P22222") == ("P11111", "P22222")


def test_duplicate_rows_and_ambiguous_source_are_not_collapsed(tmp_path):
    project = create_project(tmp_path, "integration")
    state = {
        "alpha": {"source": {"snapshot": "X", "frame": _frame(
            ("ROW000001", "1", "P1;P2", "P1;P2", "GENE", ""))}},
        "beta": {"target": {"snapshot": "Y", "frame": _frame(
            ("ROW000005", "5", "P1", "P1", "GENE", ""),
            ("ROW000009", "9", "P1", "P1", "GENE", ""),
            ("ROW000010", "10", "P2", "P2", "GENE", ""))}},
    }
    index = rebuild_cross_module_index(project, _registry(state))
    source = next(row for row in index.search("ROW000001") if row["module"] == "alpha")
    context = index.context(source["id"])
    assert context.protein_group_members == ("P1", "P2")
    matches = index.matches(context, include_other_lineages=True)
    assert {match.details["source_row"] for match in matches} == {"5", "9", "10"}
    assert all(match.status is MatchStatus.AMBIGUOUS for match in matches)


def test_rebuild_stale_and_corrupt_index(tmp_path):
    project = create_project(tmp_path, "integration")
    state = {"alpha": {"first": {"snapshot": "X", "frame": _frame(
        ("ROW000001", "1", "P1", "P1", "", ""))}}, "beta": {}}
    registry = _registry(state)
    first = open_cross_module_index(project, registry)
    assert first.healthy()
    state["alpha"]["second"] = {"snapshot": "Y", "frame": _frame(
        ("ROW000002", "2", "P2", "P2", "", ""))}
    refreshed = open_cross_module_index(project, registry)
    assert any(row["run_id"] == "second" for row in refreshed.search("P2"))
    manifest = json.loads((project.root / "analyses/integration/cross_module_index_manifest.json").read_text())
    assert manifest["schema_version"] == 1
    assert {run["snapshot"] for run in manifest["runs"]} == {"X", "Y"}
    refreshed.path.write_bytes(b"not SQLite")
    rebuilt = open_cross_module_index(project, registry)
    assert rebuilt.healthy()
    with sqlite3.connect(rebuilt.path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM records").fetchone()[0] == 2
