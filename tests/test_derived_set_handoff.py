import os
from types import SimpleNamespace

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication

from pichanalysis.core import derived_set_handoff as handoff
from pichanalysis.core.experiment_comparison import ComparisonSpec, TableSpec, compare
from pichanalysis.ui.experiment_comparison_page import ExperimentComparisonPage

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture
def comparison_project(tmp_path):
    QApplication.instance() or QApplication([])
    project = SimpleNamespace(root=tmp_path, config={"input": {"processed_file": None}})
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    pd.DataFrame({"id": ["P00001", "P00003", "P00004"]}).to_csv(a, index=False)
    pd.DataFrame({"id": ["P00002", "P00003", "P00004"]}).to_csv(b, index=False)
    compare(project, ComparisonSpec(TableSpec(str(a), "id", "uniprot"),
        TableSpec(str(b), "id", "uniprot"), "uniprot"), run_id="comparison-a")
    catalog = tmp_path / "mapping/tables/protein_catalog.csv"
    catalog.parent.mkdir(parents=True)
    pd.DataFrame({"source_row": [1, 2, 3, 4],
        "uniprot_accession": ["P00001", "P00002", "P00003", "P00004"]}).to_csv(catalog, index=False)
    root = tmp_path / "analyses/experiment_comparison/runs/comparison-a/tables"
    pd.DataFrame({"comparison_entity": ["P00004"]}).to_csv(root / "opposite_direction.csv", index=False)
    return project


@pytest.mark.parametrize("derived,destination,expected", [
    ("a_only", "go", [1]), ("b_only", "kegg", [2]),
    ("shared", "reactome", [3, 4]), ("opposite_direction", "string", [4]),
])
def test_derived_set_resolves_project_rows_without_running_analysis(
        comparison_project, derived, destination, expected):
    prepared = handoff.prepare(comparison_project, "comparison-a", derived, destination)
    assert prepared["source_rows"] == expected
    assert prepared["source"] == "Experiment Comparison"
    assert prepared["comparison_run_id"] == "comparison-a"
    assert prepared["derived_set"] == derived
    path = handoff.record(comparison_project, prepared)
    assert path.is_file()
    assert list((comparison_project.root / "analyses" / destination / "handoffs").glob("*.json"))
    assert not (comparison_project.root / "analyses" / destination / "runs").exists()


def test_unmapped_entity_is_rejected(comparison_project):
    catalog = comparison_project.root / "mapping/tables/protein_catalog.csv"
    pd.read_csv(catalog).query("uniprot_accession != 'P00001'").to_csv(catalog, index=False)
    with pytest.raises(ValueError, match="absent"):
        handoff.prepare(comparison_project, "comparison-a", "a_only", "go")


def test_gui_emits_handoff_without_running(comparison_project):
    page = ExperimentComparisonPage()
    page.set_project(comparison_project)
    page.history.setCurrentIndex(page.history.findData("comparison-a"))
    page._load_history()
    calls = []
    page.derived_target_requested.connect(calls.append)
    page.derived_set.setCurrentIndex(page.derived_set.findData("a_only"))
    page.derived_destination.setCurrentIndex(page.derived_destination.findData("go"))
    page._analyze_derived()
    assert len(calls) == 1 and calls[0]["source_rows"] == [1]
    assert not (comparison_project.root / "analyses/go").exists()
    page.close()
