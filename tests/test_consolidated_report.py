from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtCore import Qt, QSize
from PySide6.QtPdf import QPdfDocument
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication

from pichanalysis.core.consolidated_report import (
    ReportConfig, ReportError, ReportSection, SourceAdapter, discover_runs,
    generate_report, list_reports, load_report, export_report,
)
from pichanalysis.ui.consolidated_report_page import ConsolidatedReportPage


@pytest.fixture(scope="module", autouse=True)
def qt_application():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def report_fixture(tmp_path):
    project = SimpleNamespace(root=tmp_path, name="Offline project", config={})
    adapters = []
    for module in ("proteomics_qc", "differential", "go", "kegg", "string"):
        run_id = f"{module}-001"
        root = tmp_path / "analyses" / module / "runs" / run_id
        (root / "tables").mkdir(parents=True)
        (root / "plots").mkdir()
        (root / "tables" / "results.csv").write_text("gene,p_value\nA,0.02\nB,0.2\n", encoding="utf-8")
        image = QImage(40, 30, QImage.Format.Format_RGB32)
        image.fill(Qt.GlobalColor.blue)
        assert image.save(str(root / "plots" / "figure.png"))
        metadata = {"run_id": run_id, "created_at": "2026-09-24T10:00:00Z",
                    "input_lineage": module + "-input", "snapshot_id": module + "-snapshot"}
        (root / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        adapters.append(SourceAdapter(module, module, lambda _project, id=run_id: [id],
            lambda _project, _id, folder=root: {"metadata": json.loads((folder / "metadata.json").read_text()), "run_root": folder},
            lambda _project, _id, outputs: outputs["run_root"], ("tables", "plots")))
    return project, tuple(adapters)


def test_explicit_offline_report_and_history(report_fixture):
    project, adapters = report_fixture
    runs = discover_runs(project, adapters)
    assert len(runs) == 5
    config = ReportConfig(report_title="Frozen science", sections=tuple(
        ReportSection(run.module, run.run_id, tuple(a.artifact_id for a in run.artifacts), 10,
                      "User note") for run in runs), user_discussion="User discussion")
    first = generate_report(project, config, adapters)
    second = generate_report(project, config, adapters)
    assert first != second
    assert (first / "report.pdf").read_bytes().startswith(b"%PDF-")
    pdf_document = QPdfDocument()
    assert pdf_document.load(str(first / "report.pdf")) == QPdfDocument.Error.None_
    has_original_blue_figure = False
    for page_number in range(pdf_document.pageCount()):
        page_image = pdf_document.render(page_number, QSize(600, 850))
        for y in range(0, page_image.height(), 3):
            for x in range(0, page_image.width(), 3):
                color = page_image.pixelColor(x, y)
                if color.blue() > 180 and color.red() < 40 and color.green() < 40:
                    has_original_blue_figure = True
                    break
            if has_original_blue_figure:
                break
        if has_original_blue_figure:
            break
    assert has_original_blue_figure, "A selected figure must render in the PDF"
    pdf_document.close()
    assert (first / "report.html").is_file()
    manifest = json.loads((first / "report_manifest.json").read_text())
    assert manifest["status"] == "Ready"
    assert len(manifest["selected_runs"]) == 5
    assert len(manifest["selected_artifacts"]) == 10
    assert len(list_reports(project)) == 2
    assert load_report(project, first.name)["path"] == first
    html = (first / "report.html").read_text()
    assert "User discussion" in html and "User note" in html
    assert "multiple input lineages" in html
    (first / "attachments" / "section-1" / "tables" / "results.csv").write_text("tampered")
    with pytest.raises(ReportError):
        load_report(project, first.name)
    historical_html = (second / "report.html").read_bytes()
    historical_pdf = (second / "report.pdf").read_bytes()
    source_to_mutate = discover_runs(project, (adapters[0],))[0].artifacts[0].path
    source_to_mutate.write_text("changed current source", encoding="utf-8")
    assert (second / "report.html").read_bytes() == historical_html
    assert (second / "report.pdf").read_bytes() == historical_pdf
    (project.root / "current_input.csv").write_text("new current input", encoding="utf-8")
    assert load_report(project, second.name)["path"] == second
    for adapter in adapters:
        run = discover_runs(project, (adapter,))[0]
        for artifact in run.artifacts:
            artifact.path.unlink()
    assert load_report(project, second.name)["path"] == second
    exported = export_report(project, second.name, project.root / "outside.html", "html")
    assert exported.is_file()
    assert (project.root / "outside_assets").is_dir()
    assert "outside_assets/" in exported.read_text(encoding="utf-8")


def test_missing_source_blocks_generation(report_fixture):
    project, adapters = report_fixture
    run = discover_runs(project, adapters)[0]
    artifact = run.artifacts[0]
    artifact.path.unlink()
    config = ReportConfig(sections=(ReportSection(run.module, run.run_id, (artifact.artifact_id,)),))
    with pytest.raises(ReportError):
        generate_report(project, config, adapters)


def test_gui_explicit_selection(report_fixture):
    app = QApplication.instance() or QApplication([])
    project, adapters = report_fixture
    page = ConsolidatedReportPage()
    page.sources = adapters
    page.set_project(project)
    assert page.sources_tree.topLevelItemCount() == 5
    run_item = page.sources_tree.topLevelItem(0).child(0)
    run_item.setCheckState(0, Qt.CheckState.Checked)
    assert not page.generate.isEnabled()
    run_item.child(0).setCheckState(0, Qt.CheckState.Checked)
    assert page.generate.isEnabled()
    assert len(page.build_config().sections) == 1
    page._generate()
    assert page.worker is not None
    assert page.worker.wait(30000)
    app.processEvents()
    assert len(list_reports(project)) == 1


def test_two_runs_same_module_and_bounded_preview(tmp_path):
    project = SimpleNamespace(root=tmp_path, name="Two runs", config={})
    module = "differential"
    for run_id in ("A", "B"):
        folder = tmp_path / "analyses" / module / "runs" / run_id
        (folder / "tables").mkdir(parents=True)
        (folder / "tables" / "results.csv").write_text(
            "row,value\n" + "".join(f"{run_id}{index},{index}\n" for index in range(120)),
            encoding="utf-8")
        (folder / "metadata.json").write_text(json.dumps({"run_id": run_id,
            "snapshot_id": f"snapshot-{run_id}", "input_lineage": f"lineage-{run_id}"}))

    def load(_project, run_id):
        folder = tmp_path / "analyses" / module / "runs" / run_id
        return {"run_root": folder, "metadata": json.loads((folder / "metadata.json").read_text())}

    adapter = SourceAdapter(module, module, lambda _project: ["A", "B"], load,
        lambda _project, _id, outputs: outputs["run_root"], ("tables",))
    one = ReportConfig(sections=(ReportSection(module, "A", ("tables/results.csv",)),))
    first = generate_report(project, one, (adapter,))
    first_html = (first / "report.html").read_text()
    assert "A0" in first_html and "A19" in first_html and "A20" not in first_html
    assert "B0" not in first_html
    assert "Total persisted rows: 120" in first_html
    assert (first / "attachments/section-1/tables/results.csv").read_text() == (
        tmp_path / "analyses/differential/runs/A/tables/results.csv").read_text()
    both = ReportConfig(sections=(ReportSection(module, "B", ("tables/results.csv",)),
                                  ReportSection(module, "A", ("tables/results.csv",))))
    second = generate_report(project, both, (adapter,))
    manifest = json.loads((second / "report_manifest.json").read_text())
    assert [run["run_id"] for run in manifest["selected_runs"]] == ["B", "A"]
    assert manifest["selected_runs"][0]["database_snapshot"] == "snapshot-B"