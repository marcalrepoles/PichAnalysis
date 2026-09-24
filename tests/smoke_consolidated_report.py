"""Offline GUI smoke for persisted, explicitly selected report sources."""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from types import SimpleNamespace

from PySide6.QtCore import Qt
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication

from pichanalysis.core.consolidated_report import (
    SourceAdapter, export_report, list_reports, load_report,
)
from pichanalysis.ui.consolidated_report_page import ConsolidatedReportPage


def main():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory(prefix="pich_report_smoke_") as temporary:
        project = SimpleNamespace(root=Path(temporary), name="Report smoke", config={})
        adapters = []
        modules = ("proteomics_qc", "differential", "go", "kegg", "mtdna_evidence")
        for module in modules:
            run_id = module + "-A"
            folder = project.root / "analyses" / module / "runs" / run_id
            (folder / "tables").mkdir(parents=True)
            (folder / "plots").mkdir()
            (folder / "tables" / "results.csv").write_text(
                "identifier,classification,p_value,adj.P.Val\n"
                "GENE1,Not significant,NA,NA\nGENE2,Not tested,0.000001,0.02\n",
                encoding="utf-8")
            image = QImage(60, 40, QImage.Format.Format_RGB32)
            image.fill(Qt.GlobalColor.darkCyan)
            assert image.save(str(folder / "plots" / "figure.png"))
            metadata = {"run_id": run_id, "created_at": "2026-09-24T10:00:00Z",
                        "input_lineage": module + "-input", "snapshot_id": module + "-snapshot"}
            (folder / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
            adapters.append(SourceAdapter(module, module, lambda _project, id=run_id: [id],
                lambda _project, _id, path=folder: {
                    "run_root": path, "metadata": json.loads((path / "metadata.json").read_text())},
                lambda _project, _id, outputs: outputs["run_root"], ("tables", "plots")))

        page = ConsolidatedReportPage()
        page.sources = tuple(adapters)
        page.set_project(project)
        assert page.sources_tree.topLevelItemCount() == len(modules)
        for group_index in range(len(modules)):
            run_item = page.sources_tree.topLevelItem(group_index).child(0)
            run_item.setCheckState(0, Qt.CheckState.Checked)
            for artifact_index in range(run_item.childCount()):
                run_item.child(artifact_index).setCheckState(0, Qt.CheckState.Checked)
        assert page.generate.isEnabled()
        page.order.setCurrentRow(4)
        page._move(-1)
        page.notes.setPlainText("User note: α β café <script>alert(1)</script>")
        page.discussion.setPlainText("User-provided discussion only.")
        assert len(page.build_config().sections) == len(modules)
        assert "5 sections" in page.preview.text()

        page._generate()
        assert page.worker and page.worker.wait(30000)
        app.processEvents()
        first = list_reports(project)[0]["report_id"]
        report_a = load_report(project, first)
        html_a = report_a["html"].read_bytes()
        assert b"&lt;script&gt;" in html_a
        assert b"multiple input lineages" in html_a
        assert report_a["pdf"].read_bytes().startswith(b"%PDF-")

        page.title.setText("Report B")
        page._generate()
        assert page.worker and page.worker.wait(30000)
        app.processEvents()
        assert len(list_reports(project)) == 2
        assert load_report(project, first)["html"].read_bytes() == html_a

        for adapter in adapters:
            source = project.root / "analyses" / adapter.module_id / "runs" / (adapter.module_id + "-A")
            (source / "tables" / "results.csv").unlink()
            (source / "plots" / "figure.png").unlink()
        assert load_report(project, first)["html"].read_bytes() == html_a
        assert report_a["pdf"].read_bytes().startswith(b"%PDF-")
        export = project.root / "exports"
        export.mkdir()
        assert export_report(project, first, export / "report.html", "html").is_file()
        assert (export / "report_assets").is_dir()
        assert export_report(project, first, export / "report.pdf", "pdf").is_file()
        assert export_report(project, first, export / "bundle", "bundle").is_dir()
        assert load_report(project, first)["manifest"]["selected_runs"][0]["run_id"]
        print("Consolidated Report offline GUI smoke: PASS")


if __name__ == "__main__":
    main()
