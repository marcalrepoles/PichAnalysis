"""Offline, frozen-run safeguards for Biological Context."""
from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest
from PySide6.QtGui import QImage

from pichanalysis.core import biological_context as bc
from pichanalysis.core.entity_identity import EntityContext, InputLineage, MatchStatus


def context(module="differential", run="A", row=1, accession="P12345-2"):
    return EntityContext(module, run, "all_results", "F1", row, accession,
        uniprot_accessions=(accession,), gene_symbols=("GENE1",),
        input_lineage=InputLineage("project", "hash-A"))


def test_invalid_and_ambiguous_source(monkeypatch):
    class Index:
        def run_status(self, module, run): return ("Ready", "")
        def search(self, query):
            return [{"module": "differential", "run_id": "A", "record_type": "all_results",
                     "source_row": row, "id": row} for row in (1, 2)]
        def context(self, row): return context(row=row)
    monkeypatch.setattr(bc, "open_cross_module_index", lambda project: Index())
    with pytest.raises(bc.BiologicalContextError, match="valid"):
        bc.resolve_entity(None, "differential", "A", "")
    with pytest.raises(bc.BiologicalContextError, match="specific source row"):
        bc.resolve_entity(None, "differential", "A", "F1")
    _, result = bc.resolve_entity(None, "differential", "A", "F1", 2)
    assert result.source_row == 2


def test_incomplete_source_cannot_borrow_other_run(monkeypatch):
    index = SimpleNamespace(run_status=lambda *args: ("Incomplete", "missing"))
    monkeypatch.setattr(bc, "open_cross_module_index", lambda project: index)
    with pytest.raises(bc.BiologicalContextError, match="incomplete"):
        bc.resolve_entity(None, "differential", "A", "F1")


def test_isoform_and_source_row_are_preserved():
    entity = context(row=9, accession="P12345-2")
    assert entity.source_row == 9
    assert entity.uniprot_accessions == ("P12345-2",)
    assert entity.uniprot_accessions != ("P12345",)


@pytest.mark.parametrize(("classification", "direction", "expected"), [
    ("TRUE", "higher_in_condition_A", "up"),
    ("TRUE", "higher_in_condition_B", "down"),
    ("FALSE", "higher_in_condition_A", "none"),
    ("TRUE", "no_difference", "none"),
])
def test_only_persisted_significant_direction(monkeypatch, classification, direction, expected):
    frame = pd.DataFrame([{"source_row": "1", "feature_id": "F1", "uniprot_accession": "P12345-2",
        "gene_symbol": "GENE1", "combined_significant": classification,
        "effect_direction": direction}])
    monkeypatch.setattr(bc.differential_analysis, "load_run",
        lambda project, run: {"tables": {"all_results": frame}})
    assert bc.frozen_direction(None, "A", context(), {"UniProt": "P12345-2"}) == expected
    assert bc.frozen_direction(None, None, context(), {"UniProt": "P12345-2"}) == "none"
    assert bc.frozen_direction(None, "A", context(), {"UniProt": "OTHER"}) == "none"


def test_explicit_differential_choices_never_take_latest(tmp_path):
    import sqlite3
    database = tmp_path / "index.sqlite"
    with sqlite3.connect(database) as db:
        db.execute("CREATE TABLE runs(module TEXT,run_id TEXT,status TEXT,input_hash TEXT)")
        db.executemany("INSERT INTO runs VALUES (?,?,?,?)", [
            ("differential", "A", "Ready", "hash-A"),
            ("differential", "B", "Ready", "hash-A"),
            ("differential", "C", "Ready", "other"),
            ("differential", "D", "Incomplete", "hash-A")])
    choices = bc.differential_choices(SimpleNamespace(path=database), context(module="mapping"))
    assert choices == ["A", "B"]
    assert bc.other_differential_runs(SimpleNamespace(path=database), choices) == ["C"]


def test_official_images_are_snapshot_specific(tmp_path):
    root = tmp_path / "databases"
    old = root / "kegg/hsa/snapshots/X/images"
    new = root / "kegg/hsa/snapshots/Y/images"
    old.mkdir(parents=True); new.mkdir(parents=True)
    (old / "hsa00010.png").write_bytes(b"old")
    (new / "hsa00010.png").write_bytes(b"new")
    assert bc._official("kegg", "X", "hsa00010", root) == old / "hsa00010.png"
    assert bc._official("kegg", "missing", "hsa00010", root) is None


def test_reactome_official_diagram_is_local_and_missing_is_safe(tmp_path):
    root = tmp_path / "databases"
    snapshot = root / "reactome/human/snapshots/X"
    (snapshot / "diagrams").mkdir(parents=True)
    (snapshot / "diagram_index.csv").write_text(
        "Reactome_ID,relative_path\nR-HSA-1,R-HSA-1.png\n", encoding="utf-8")
    (snapshot / "diagrams/R-HSA-1.png").write_bytes(b"diagram")
    assert bc._official("reactome", "X", "R-HSA-1", root) == snapshot / "diagrams/R-HSA-1.png"
    assert bc._official("reactome", "Y", "R-HSA-1", root) is None


def test_string_reads_only_persisted_direct_edges(monkeypatch):
    nodes = pd.DataFrame({"string_protein_id": ["S1", "S2", "S3"], "gene_symbol": ["A", "B", "C"]})
    edges = pd.DataFrame({"protein_a": ["S1", "S2"], "protein_b": ["S2", "S3"], "combined_score": [900, 800]})
    output = {"tables": {"expanded_nodes": nodes, "expanded_edges": edges},
              "metadata": {"snapshot_id": "X", "network_type": "physical"}}
    match = SimpleNamespace(target_module="string", target_run_id="A", status=MatchStatus.MAPPED,
        details={"record_id": "1"})
    index = SimpleNamespace(matches=lambda *args, **kwargs: [match],
        run_status=lambda *args: ("Ready", ""), path=None)
    monkeypatch.setattr(bc, "_details", lambda *args: {"string_protein_id": "S1"})
    monkeypatch.setattr(bc.string_analysis, "read_string_outputs", lambda *args: output)
    items = bc.find_contexts(None, index, context(module="mapping"))
    assert len(items) == 1
    assert set(items[0].members.string_protein_id) == {"S1", "S2"}
    assert len(items[0].edges) == 1


def test_no_string_network_is_not_invented(monkeypatch):
    index = SimpleNamespace(matches=lambda *args, **kwargs: [], path=None)
    assert bc.find_contexts(None, index, context(module="mapping")) == []


def test_generated_history_retains_original_run_and_snapshot(tmp_path):
    project = SimpleNamespace(root=tmp_path)
    item = bc.ContextItem("kegg", "A", "hsa00010", "Pathway", "X", "hsa:1",
        pd.DataFrame({"kegg_gene_id": ["hsa:1"]}), pd.DataFrame())
    image = QImage(10, 10, QImage.Format_RGB32)
    image.fill(0xFFFFFF)
    root = bc.save_context_figure(project, context(), item, "A", image)
    metadata, figure = bc.load_saved_context(root)
    assert figure.is_file()
    assert (metadata["source_run"], metadata["item_run"], metadata["snapshot"]) == ("A", "A", "X")
    assert metadata["legend"] == {"up": bc.GREEN, "down": bc.RED, "no_direction": bc.GRAY}
    assert bc.list_saved_contexts(project) == [root]

def test_render_selected_entity_and_export_png(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication
    from pichanalysis.core.biological_context_render import render_context
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    application = QApplication.instance() or QApplication([])
    project = SimpleNamespace(root=tmp_path)
    item = bc.ContextItem("string", "A", "S1", "Network", "X", "S1",
        pd.DataFrame({"string_protein_id": ["S1", "S2"], "gene_symbol": ["A", "B"]}),
        pd.DataFrame({"protein_a": ["S1"], "protein_b": ["S2"]}))
    image = render_context(project, context(), item)
    assert image.width() >= 1000 and not image.isNull()
    output = tmp_path / "context.png"
    assert image.save(str(output), "PNG") and output.is_file()

def test_differential_kegg_reactome_runs_and_snapshots_stay_separate(tmp_path, monkeypatch):
    project = SimpleNamespace(root=tmp_path)
    for run, identifier in (("A", "hsa00010"), ("B", "hsa00020")):
        folder = tmp_path / run / "mapping"
        folder.mkdir(parents=True)
        pd.DataFrame({"pathway_id": [identifier], "kegg_gene_id": ["hsa:1"],
                      "name": [f"Pathway {run}"]}).to_csv(folder / "pathway_membership.csv", index=False)
    matches = [SimpleNamespace(target_module=module, target_run_id=run,
        status=MatchStatus.MAPPED, details={"record_id": "1"})
        for module, run in (("kegg", "A"), ("kegg", "B"), ("reactome", "R"))]
    index = SimpleNamespace(matches=lambda *args, **kwargs: matches,
        run_status=lambda *args: ("Ready", ""), path=None)
    monkeypatch.setattr(bc, "_details", lambda index, project, match:
        {"kegg_gene_id": "hsa:1"} if match.target_module == "kegg" else
        {"reactome_entity_key": "U:P12345-2"})
    monkeypatch.setattr(bc.kegg_analysis, "read_kegg_outputs", lambda project, run:
        SimpleNamespace(root=tmp_path / run, metadata={"snapshot_id": f"snap-{run}"}))
    membership = pd.DataFrame({"Reactome_ID": ["R-HSA-1"], "Pathway": ["Reactome fixture"],
        "reactome_entity_key": ["U:P12345-2"], "UniProt": ["P12345-2"]})
    monkeypatch.setattr(bc.reactome_analysis, "read_reactome_outputs", lambda project, run:
        SimpleNamespace(membership=membership, metadata={"snapshot_id": "snap-R"}))
    monkeypatch.setattr(bc, "_official", lambda module, snapshot, item_id, database_root=None: None)
    found = bc.find_contexts(project, index, context(), include_other_lineages=True)
    assert {(item.module, item.run_id, item.item_id, item.snapshot) for item in found} == {
        ("kegg", "A", "hsa00010", "snap-A"),
        ("kegg", "B", "hsa00020", "snap-B"),
        ("reactome", "R", "R-HSA-1", "snap-R")}


def test_multiple_target_identities_are_not_silently_chosen(monkeypatch):
    matches = [SimpleNamespace(target_module="string", target_run_id="S1",
        status=MatchStatus.MAPPED, details={"record_id": str(index)}) for index in (1, 2)]
    index = SimpleNamespace(matches=lambda *args, **kwargs: matches,
        run_status=lambda *args: ("Ready", ""), path=None)
    monkeypatch.setattr(bc, "_details", lambda index, project, match:
        {"string_protein_id": "S" + match.details["record_id"]})
    assert bc.find_contexts(None, index, context()) == []

def test_module_specific_identifier_fallback(monkeypatch):
    class Index:
        def run_status(self, module, run): return ("Ready", "")
        def search(self, query):
            return ([{"module": "kegg", "run_id": "K1", "record_type": "mapping",
                      "source_row": 1, "id": 7}] if query == "kegg_gene_id:hsa:1" else [])
        def context(self, record_id): return context(module="kegg", run="K1")
    monkeypatch.setattr(bc, "open_cross_module_index", lambda project: Index())
    _, result = bc.resolve_entity(None, "kegg", "K1", "hsa:1")
    assert result.source_run_id == "K1"
