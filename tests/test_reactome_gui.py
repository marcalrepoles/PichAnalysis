import base64
import json
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication, QFileDialog
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QMessageBox

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.reactome_analysis import (
    MissingReactomeOutputError, ReactomeParameters, ReactomeRExecutionError, TargetOutsideBackgroundError,
    list_reactome_runs, reactome_pathway_proteins, reactome_protein_pathways,
    read_reactome_outputs,
)
from pichanalysis.ui.main_window import MainWindow
from pichanalysis.ui.reactome_page import ReactomeAnalysisWorker, ReactomePage


PNG = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def project_fixture(tmp_path, human=True, presence=True):
    project = create_project(tmp_path, "GuiProject")
    set_organism(project, "Homo sapiens" if human else "Mus musculus", "9606" if human else "10090")
    catalog = pd.DataFrame([
        {"source_row": 1, "original_id": "P1", "mapping_status": "mapped_unique", "uniprot_accession": "U1", "ncbi_gene_id": "N1", "gene_symbol": "G1", "protein_name": "Protein one"},
        {"source_row": 2, "original_id": "P2", "mapping_status": "mapped_unique", "uniprot_accession": "U2", "ncbi_gene_id": "N2", "gene_symbol": "G2", "protein_name": "Protein two"},
    ])
    catalog.to_csv(project.root / "mapping" / "tables" / "protein_catalog.csv", index=False)
    if presence:
        directory = project.root / "analyses" / "presence_absence" / "tables"
        directory.mkdir(parents=True)
        pd.DataFrame({"source_row": [1, 2], "classification": ["A-specific", "Shared"]}).to_csv(directory / "classification.csv", index=False)
    return project


def reactome_source(root):
    root.mkdir()
    values = {
        "ReactomePathways.txt": "R-HSA-1\tPathway one\tHomo sapiens\n",
        "ReactomePathwaysRelation.txt": "R-HSA-1\tR-HSA-2\n",
        "UniProt2Reactome.txt": "U1\tR-HSA-1\thttps://x\tPathway one\tIEA\tHomo sapiens\n",
        "NCBI2Reactome.txt": "N1\tR-HSA-1\thttps://x\tPathway one\tIEA\tHomo sapiens\n",
        "humanPathwaysWithDiagrams.txt": "R-HSA-1\n",
        "pathway2summation.txt": "R-HSA-1\tSummary\n",
    }
    for name, value in values.items():
        (root / name).write_text(value, encoding="utf-8")
    return root


def ready_manager(tmp_path):
    manager = DatabaseManager(tmp_path / "databases")
    snapshot = manager.reactome.install_from_directory(reactome_source(tmp_path / "source"))
    manifest = manager.reactome.manifest(snapshot)
    manifest["release_version"] = "97"
    manager.reactome._write_manifest(snapshot, manifest)
    return manager


def write_run(project, run_id, marker, graphs=True, mandatory=True):
    root = project.root / "analyses" / "Reactome" / "runs" / run_id
    for folder in ("mapping", "pathways", "frequency", "enrichment", "graphs"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    mapping = pd.DataFrame([
        {"source_row": 1, "original_id": f"{marker}-P1", "reactome_entity_key": "UP:U1", "uniprot_accession": "U1", "ncbi_gene_id": "N1", "gene_symbol": f"{marker}-G1", "protein_name": "Protein one", "mapping_route": "UniProt", "mapping_status": "mapped_unique"},
        {"source_row": 2, "original_id": f"{marker}-P2", "reactome_entity_key": "UP:U2", "uniprot_accession": "U2", "ncbi_gene_id": "N2", "gene_symbol": f"{marker}-G2", "protein_name": "Protein two", "mapping_route": "UniProt", "mapping_status": "mapped_unique"},
    ])
    mapping.to_csv(root / "mapping" / "reactome_entity_mapping.csv", index=False)
    mapping.iloc[0:0].to_csv(root / "mapping" / "unmapped.csv", index=False)
    mapping.iloc[0:0].to_csv(root / "mapping" / "ambiguous.csv", index=False)
    membership = pd.DataFrame([
        {"Reactome_ID": "R-HSA-1", "Pathway": f"{marker} Pathway one", "reactome_entity_key": "UP:U1", "Gene_symbol": f"{marker}-G1", "UniProt": "U1", "NCBI_Gene_ID": "N1", "Original_ID": f"{marker}-P1", "Mapping_route": "UniProt"},
        {"Reactome_ID": "R-HSA-2", "Pathway": f"{marker} Pathway two", "reactome_entity_key": "UP:U2", "Gene_symbol": f"{marker}-G2", "UniProt": "U2", "NCBI_Gene_ID": "N2", "Original_ID": f"{marker}-P2", "Mapping_route": "UniProt"},
    ])
    membership.to_csv(root / "pathways" / "pathway_membership.csv", index=False)
    pd.DataFrame({"parent_pathway_id": ["R-HSA-1"], "parent_pathway_name": [f"{marker} Parent"], "child_pathway_id": ["R-HSA-2"], "child_pathway_name": [f"{marker} Child"]}).to_csv(root / "pathways" / "pathway_hierarchy.csv", index=False)
    pd.DataFrame({"pathway_id": ["R-HSA-2"], "top_level_ancestors": ["R-HSA-1"], "minimum_depth": [1], "cycle_detected": [False]}).to_csv(root / "pathways" / "pathway_ancestry.csv", index=False)
    pd.DataFrame({"Reactome_ID": ["R-HSA-1"], "Pathway": [f"{marker} Pathway one"], "Protein_count": [1], "Protein_fraction": [0.5]}).to_csv(root / "frequency" / "reactome_frequency.csv", index=False)
    enrichment = pd.DataFrame({"Reactome_ID": ["R-HSA-1"], "Pathway": [f"{marker} Pathway one"], "Target_count": [1], "Target_size": [1], "Background_count": [1], "Background_size": [2], "GeneRatio": [1.0], "BgRatio": [0.5], "p_value": [0.5], "FDR": [0.5], "Entities": ["UP:U1"]})
    enrichment.to_csv(root / "enrichment" / "reactome_enrichment_all.csv", index=False)
    enrichment.iloc[0:0].to_csv(root / "enrichment" / "reactome_enrichment_significant.csv", index=False)
    enrichment.iloc[0:0].to_csv(root / "enrichment" / "enrichment_excluded.csv", index=False)
    pd.DataFrame({"metric": ["Snapshot ID", "Input entities"], "value": [f"{marker}-snapshot", 2]}).to_csv(root / "summary.csv", index=False)
    metadata = {"run_id": run_id, "target_definition": f"{marker} target", "background_definition": f"{marker} background", "snapshot_id": f"{marker}-snapshot", "reactome_release": marker}
    (root / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    (root / "Reactome_analysis.xlsx").write_bytes(f"{marker}-workbook".encode())
    if graphs:
        for stem in ("reactome_frequency", "reactome_enrichment_dot", "reactome_enrichment_fdr"):
            (root / "graphs" / f"{stem}.png").write_bytes(PNG)
            (root / "graphs" / f"{stem}.pdf").write_bytes(b"%PDF-1.4 fixture")
    if not mandatory:
        (root / "frequency" / "reactome_frequency.csv").unlink()
    return root


def test_loader_lists_runs_and_optional_graphs(tmp_path):
    project = project_fixture(tmp_path)
    write_run(project, "B", "beta")
    write_run(project, "A", "alpha", graphs=False)
    assert list_reactome_runs(project) == ["A", "B"]
    output = read_reactome_outputs(project, "A")
    assert output.graphs == () and output.summary.iloc[0]["value"] == "alpha-snapshot"


def test_loader_missing_mandatory_artifact(tmp_path):
    project = project_fixture(tmp_path)
    write_run(project, "broken", "broken", mandatory=False)
    with pytest.raises(MissingReactomeOutputError, match="reactome_frequency.csv"):
        read_reactome_outputs(project, "broken")


def test_navigation_helpers_are_local_and_keyed_by_ids(tmp_path):
    project = project_fixture(tmp_path)
    write_run(project, "A", "alpha")
    outputs = read_reactome_outputs(project, "A")
    proteins = reactome_pathway_proteins(outputs, "R-HSA-1")
    assert list(proteins.reactome_entity_key) == ["UP:U1"]
    assert proteins.iloc[0]["protein_name"] == "Protein one"
    pathways = reactome_protein_pathways(outputs, "UP:U2")
    assert list(pathways.Reactome_ID) == ["R-HSA-2"]
    assert reactome_pathway_proteins(outputs, "R-HSA-X").empty
    assert reactome_protein_pathways(outputs, "UP:MISSING").empty


def test_page_readiness_targets_background_and_advanced(app, tmp_path):
    project = project_fixture(tmp_path, presence=False)
    page = ReactomePage(DatabaseManager(tmp_path / "empty"))
    page.set_project(project)
    assert "not installed" in page.readiness.text().lower()
    assert not page.run_button.isEnabled()
    assert page.target.findText("Shared") >= 0
    assert not page.target.model().item(page.target.findText("Shared")).isEnabled()
    assert page.background.currentData() == "All mapped entities"
    assert page.fdr.value() == pytest.approx(0.05)
    assert page.minimum.value() == 3 and page.top_n.value() == 20

    manager = ready_manager(tmp_path)
    page.manager = manager
    directory = project.root / "analyses" / "presence_absence" / "tables"
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"source_row": [1, 2], "classification": ["A-specific", "Shared"]}).to_csv(directory / "classification.csv", index=False)
    page.set_project(project)
    assert page.run_button.isEnabled()
    assert page.target.findText("A-specific") >= 0
    assert page.target.model().item(page.target.findText("Shared")).isEnabled()
    page.close()


def test_page_unsupported_organism(app, tmp_path):
    project = project_fixture(tmp_path, human=False)
    page = ReactomePage(ready_manager(tmp_path))
    page.set_project(project)
    assert "Homo sapiens only" in page.readiness.text()
    assert not page.run_button.isEnabled()
    page.close()


def test_page_results_navigation_graphs_and_exports(app, tmp_path, monkeypatch):
    project = project_fixture(tmp_path)
    manager = ready_manager(tmp_path)
    write_run(project, "A", "alpha")
    outputs = read_reactome_outputs(project, "A")
    page = ReactomePage(manager)
    page.set_project(project)
    page.show_outputs(outputs)
    assert page.frequency.rowCount() == 1 and page.enrichment.rowCount() == 1
    assert page.no_significant.isVisible() is False or page.no_significant.text().startswith("No pathways")
    assert page.pathway_members.rowCount() == 1
    assert page.protein_pathways.rowCount() == 1
    assert page.graph_choice.count() == 3 and page.graph_preview.pixmap() is not None

    table_destination = tmp_path / "frequency-export.csv"
    page.tabs.setCurrentIndex(2)
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args: (str(table_destination), "CSV (*.csv)"))
    page._export_current_table()
    assert pd.read_csv(table_destination).iloc[0]["Reactome_ID"] == "R-HSA-1"
    workbook_destination = tmp_path / "export.xlsx"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args: (str(workbook_destination), "Excel workbook (*.xlsx)"))
    page._export_workbook()
    assert workbook_destination.read_bytes() == outputs.workbook.read_bytes()
    graph_destination = tmp_path / "plot.pdf"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args: (str(graph_destination), "PDF (*.pdf)"))
    page._export_graph()
    assert graph_destination.read_bytes().startswith(b"%PDF")
    page.close()


def test_history_a_b_a_and_incomplete_run(app, tmp_path):
    project = project_fixture(tmp_path)
    manager = ready_manager(tmp_path)
    write_run(project, "A", "alpha")
    write_run(project, "B", "beta")
    write_run(project, "broken", "broken", mandatory=False)
    page = ReactomePage(manager)
    page.set_project(project)
    for run_id, marker in (("A", "alpha"), ("B", "beta"), ("A", "alpha")):
        page.history.setCurrentIndex(page.history.findData(run_id))
        app.processEvents()
        assert page.outputs.metadata["snapshot_id"] == f"{marker}-snapshot"
        assert marker in str(page.summary.item(0, 1).text())
        assert marker in page.frequency.item(0, 1).text()
        assert marker in page.enrichment.item(0, 1).text()
        assert marker in page.mapping.item(0, 1).text()
        assert marker in page.pathway_members.item(0, 1).text()
    page.history.setCurrentIndex(page.history.findData("broken"))
    app.processEvents()
    assert page.outputs is None and "incomplete" in page.history_status.text().lower()
    page.close()


def test_worker_success_target_error_and_general_error(app, tmp_path):
    project = project_fixture(tmp_path)
    manager = ready_manager(tmp_path)
    write_run(project, "A", "alpha")
    outputs = read_reactome_outputs(project, "A")
    success = []
    worker = ReactomeAnalysisWorker(project, manager, object(), "x", ReactomeParameters(), runner=lambda *a, **k: outputs)
    worker.succeeded.connect(success.append)
    worker.run()
    assert success == [outputs]

    outside = []
    def outside_runner(*args, **kwargs):
        raise TargetOutsideBackgroundError({"entities_outside_background": ["UP:U2"]})
    worker = ReactomeAnalysisWorker(project, manager, object(), "y", ReactomeParameters(), runner=outside_runner)
    worker.target_outside_background.connect(outside.append)
    worker.run()
    assert outside[0]["entities_outside_background"] == ["UP:U2"]

    failures = []
    worker = ReactomeAnalysisWorker(project, manager, object(), "z", ReactomeParameters(), runner=lambda *a, **k: (_ for _ in ()).throw(MissingReactomeOutputError("Missing expected output")))
    worker.failed.connect(failures.append)
    worker.run()
    assert failures == ["Missing expected output"]
    failures.clear()
    worker = ReactomeAnalysisWorker(project, manager, object(), "r", ReactomeParameters(), runner=lambda *a, **k: (_ for _ in ()).throw(ReactomeRExecutionError("R execution failure: planned")))
    worker.failed.connect(failures.append)
    worker.run()
    assert failures == ["R execution failure: planned"]


def test_main_window_target_background_cancel_and_continue(app, tmp_path, monkeypatch):
    project = project_fixture(tmp_path)
    window = MainWindow()
    window.project = project
    window.analyses_page.reactome_page.set_project(project)
    calls = []
    monkeypatch.setattr(window, "_run_reactome", lambda parameters, allow=False: calls.append((parameters, allow)))
    parameters = {"target_selection": "Shared"}
    monkeypatch.setattr(window, "_confirm_reactome_adjustment", lambda details: False)
    window.reactome_worker = object()
    window._reactome_target_outside({"entities_outside_background": ["UP:U2"]}, parameters)
    assert calls == []
    monkeypatch.setattr(window, "_confirm_reactome_adjustment", lambda details: True)
    window.reactome_worker = object()
    window._reactome_target_outside({"entities_outside_background": ["UP:U2"]}, parameters)
    assert calls == [(parameters, True)]
    window.close()


def test_python_r_gui_end_to_end_smoke(app, tmp_path):
    project = project_fixture(tmp_path)
    manager = ready_manager(tmp_path)
    page = ReactomePage(manager)
    page.set_project(project)
    page.minimum.setValue(1)
    page.fdr.setValue(0.05)
    page.top_n.setValue(20)
    assert page.run_button.isEnabled()
    parameters = ReactomeParameters(**page.parameters())
    runtime = RRuntime()
    assert runtime.available
    successes, failures = [], []
    worker = ReactomeAnalysisWorker(project, manager, runtime, "gui-smoke", parameters)
    worker.succeeded.connect(lambda outputs: (successes.append(outputs), page.show_outputs(outputs)))
    worker.failed.connect(failures.append)
    page.set_running(True)
    worker.run()
    page.set_running(False)
    assert not failures and len(successes) == 1
    assert page.outputs.metadata["run_id"] == "gui-smoke"
    assert page.frequency.rowCount() >= 1
    assert page.enrichment.rowCount() >= 1
    assert page.pathway_members.rowCount() >= 1
    assert page.protein_pathways.rowCount() >= 1
    assert page.graph_choice.count() == 3
    assert page.history.findData("gui-smoke") >= 0
    assert (project.root / "analyses" / "Reactome" / "runs" / "gui-smoke" / "Reactome_analysis.xlsx").is_file()
    page.close()


def test_application_close_is_blocked_while_reactome_worker_runs(app, monkeypatch):
    window = MainWindow()

    class RunningWorker:
        @staticmethod
        def isRunning():
            return True

    window.reactome_worker = RunningWorker()
    monkeypatch.setattr(QMessageBox, "information", lambda *args: QMessageBox.StandardButton.Ok)
    event = QCloseEvent()
    window.closeEvent(event)
    assert not event.isAccepted()
    window.reactome_worker = None
    window.close()
