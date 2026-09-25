import json
import os
from types import SimpleNamespace

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication

from pichanalysis.core.experiment_comparison import ComparisonSpec, TableSpec, compare as compare_base
from pichanalysis.core import functional_comparison as fc


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def side(ids, *, names=None, counts=None, members=None, significant=()):
    names = names or ids
    counts = counts or [1] * len(ids)
    members = members or ["P12345"] * len(ids)
    frequency = pd.DataFrame({"pathway_id": ids, "name": names, "gene_count": counts,
        "gene_fraction": ["0.5"] * len(ids), "genes": members})
    enrichment = pd.DataFrame({"pathway_id": ids, "p_value": ["0.01"] * len(ids),
        "p_adjust": ["0.02"] * len(ids), "target_gene_count": counts})
    return {"frequency": frequency, "enrichment": enrichment,
        "significant": pd.DataFrame({"pathway_id": list(significant)})}


def test_kegg_stable_id_name_change_members_and_enrichment():
    a = side(["hsa00010", "hsa00020"], names=["Old name", "A"],
        members=["P12345;P12346", "P99999"], significant=["hsa00010"])
    b = side(["hsa00010", "hsa00030"], names=["New name", "B"],
        members=["P12346;Q11111", "B22222"], significant=["hsa00010"])
    master, members = fc._compare_unit("kegg", "Pathways", a, b, "snapshot-a", "snapshot-b")
    assert set(master.entity_id) == {"hsa00010", "hsa00020", "hsa00030"}
    shared = master.set_index("entity_id").loc["hsa00010"]
    assert shared.name_A == "Old name" and shared.name_B == "New name"
    assert shared.comparison_class == "Enriched in both"
    assert shared.A_FDR == "0.02" and shared.B_FDR == "0.02"
    detail = members[members.entity_id == "hsa00010"].set_index("member_id")
    assert detail.loc["P12346", "member_class"] == "Shared"
    assert detail.loc["P12345", "member_class"] == "A only"
    assert detail.loc["Q11111", "member_class"] == "B only"


def test_go_ontology_ids_and_evidence_are_not_combined():
    a = {"frequency": pd.DataFrame({"GO_ID": ["GO:0001"], "Description": ["alpha"],
        "Protein_count": [1]}), "membership": pd.DataFrame({"GO_ID": ["GO:0001"],
        "entity_id": ["GENE1"]})}
    b = {"frequency": pd.DataFrame({"GO_ID": ["GO:0001"], "Description": ["beta"],
        "Protein_count": [1]}), "membership": pd.DataFrame({"GO_ID": ["GO:0001"],
        "entity_id": ["GENE2"]})}
    master, members = fc._compare_unit("go", "BP", a, b, "", "")
    assert len(master) == 1 and master.iloc[0].name_A == "alpha" and master.iloc[0].name_B == "beta"
    assert set(members.member_class) == {"A only", "B only"}
    assert fc.UNITS["go"] == ("BP", "MF", "CC")


def test_reactome_stable_id_and_hierarchy():
    a = {"frequency": pd.DataFrame({"Reactome_ID": ["R-HSA-1"], "Pathway": ["old"],
        "Protein_count": [1]}), "hierarchy": pd.DataFrame({"Reactome_ID": ["R-HSA-1"],
        "parent_pathway": ["R-HSA-0"]})}
    b = {"frequency": pd.DataFrame({"Reactome_ID": ["R-HSA-1"], "Pathway": ["new"],
        "Protein_count": [1]})}
    master, _ = fc._compare_unit("reactome", "Pathways", a, b, "a", "b")
    assert "R-HSA-0" in master.iloc[0].A_details
    assert master.iloc[0].snapshot_A == "a" and master.iloc[0].snapshot_B == "b"


@pytest.mark.parametrize("module,unit,identifier,column", [
    ("mitocarta", "Subcompartments", "Matrix", "Subcompartment"),
    ("mitocarta", "MitoPathways", "ATP", "MitoPathway"),
    ("mtdna", "Categories", "mtDNA encoded", "category"),
    ("mtdna", "Sources", "MitoCarta", "source"),
    ("domains", "InterPro", "IPR0001", "InterPro_ID"),
    ("domains", "Pfam", "PF0001", "Pfam_ID"),
    ("domains", "Pfam architectures", "PF1 > PF2", "Architecture"),
    ("complexes", "Complexes", "CPX-1", "complex_id"),
])
def test_module_units_are_separate(module, unit, identifier, column):
    frame = pd.DataFrame({column: [identifier], "Gene_count": [1],
        "Protein_count": [1], "target_member_count": [1], "target_entity_count": [1],
        "Genes": ["GENE1"], "Proteins": ["P12345"], "target_member_proteins": ["P12345"]})
    master, _ = fc._compare_unit(module, unit, {"frequency": frame}, {"frequency": frame}, "", "")
    assert len(master) == 1
    assert master.iloc[0].entity_id == identifier


def test_string_undirected_edges_and_nodes():
    edges_a = pd.DataFrame({"protein_a": ["9606.A", "9606.C"],
        "protein_b": ["9606.B", "9606.D"], "combined_score": [900, 800]})
    edges_b = pd.DataFrame({"protein_a": ["9606.B", "9606.E"],
        "protein_b": ["9606.A", "9606.F"], "combined_score": [850, 750]})
    master, _ = fc._compare_unit("string", "Edges", {"frequency": edges_a},
        {"frequency": edges_b}, "a", "b")
    assert set(master.entity_id) == {"9606.A|9606.B", "9606.C|9606.D", "9606.E|9606.F"}
    assert master.set_index("entity_id").loc["9606.A|9606.B", "presence_class"] == "Shared"
    assert fc.UNITS["string"] == ("Nodes", "Edges")


def test_persisted_history_and_snapshot_warning(tmp_path, monkeypatch):
    app = QApplication.instance() or QApplication([])
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    for path in (a, b):
        pd.DataFrame({"id": ["P12345"], "description": ["protein"]}).to_csv(path, index=False)
    project = SimpleNamespace(root=tmp_path)
    compare_base(project, ComparisonSpec(TableSpec(str(a), "id", "uniprot"),
        TableSpec(str(b), "id", "uniprot"), "uniprot"), run_id="comparison-a")
    sources = {"run-a": ({"snapshot_id": "old", "minimum_overlap": 2},
        {"Pathways": side(["hsa00010"], significant=["hsa00010"])}),
        "run-b": ({"snapshot_id": "new", "minimum_overlap": 3},
        {"Pathways": side(["hsa00010"], significant=[])})}
    monkeypatch.setattr(fc, "_source", lambda project, module, run_id: sources[run_id])
    result_root = fc.compare(project, fc.FunctionalSpec("comparison-a", "kegg", "run-a", "run-b"))
    loaded = fc.load_functional(project, "comparison-a", "kegg", result_root.name)
    assert "different database snapshots" in " ".join(loaded["config"]["warnings"])
    assert "minimum_overlap" in " ".join(loaded["config"]["warnings"])
    assert loaded["tables"]["Pathways"]["master"].iloc[0].comparison_class == "Enriched only in A"
    assert (result_root / "functional_comparison.xlsx").is_file()
    sources.clear()
    assert len(fc.load_functional(project, "comparison-a", "kegg", result_root.name)["tables"]["Pathways"]["master"]) == 1
