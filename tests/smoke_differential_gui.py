"""Offscreen GUI -> Python worker -> R -> frozen history smoke."""
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication

from pichanalysis.ui.differential_analysis_page import DifferentialAnalysisPage
from smoke_differential_preparation import fixture_project


def wait_for_run(page):
    loop = QEventLoop()
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(loop.quit)
    page.worker.finished.connect(loop.quit)
    timer.start(180000)
    loop.exec()
    assert timer.isActive(), "Differential worker timed out"
    QApplication.processEvents()
    assert not page.is_running(), page.status.text()
    assert page.status.text() not in ("Preparing differential-analysis data...", "Running differential statistics...")


def smoke():
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory() as folder:
        project, columns, source = fixture_project(Path(folder))
        data = pd.read_csv(source)
        data.loc[data.ProteinGroup == "P11", ["B1", "B2", "B3"]] = [20, 21, 22]
        data.loc[data.ProteinGroup == "P12", ["B1", "B2", "B3"]] = [55, 56, 57]
        data["Gene"] = [f"GENE{i:02d}" for i in range(len(data))]
        data["UniProt"] = [f"P{i:05d}" for i in range(len(data))]
        project.config["columns"]["Gene"] = {"role": "identifier", "identifier_type": "gene_symbol"}
        project.config["columns"]["UniProt"] = {"role": "identifier", "identifier_type": "uniprot"}
        project.save()
        data.to_csv(source, index=False)
        page = DifferentialAnalysisPage()
        page.set_project(project)
        assert page.stage_tabs.count() == 3
        assert page.condition_a.currentText() == "A" and page.condition_b.currentText() == "B"
        assert page.samples.rowCount() == 6
        assert page.prepare_button.isEnabled()
        assert page.imputation.currentData() == "none"
        page._start_preparation(); wait_for_run(page)
        assert page.preparation is not None, page.status.text()
        prep_a = page.preparation["metadata"]["run_id"]
        eligibility = page.preparation["tables"]["feature_eligibility"]
        assert eligibility.loc[4, "pattern"] == "detected_only_condition_A"
        assert eligibility.loc[5, "pattern"] == "detected_only_condition_B"
        assert eligibility.loc[6, "pattern"] == "sparse_only_condition_A"
        assert eligibility.loc[7, "pattern"] == "all_missing"
        assert eligibility.loc[0, "feature_id"] != eligibility.loc[1, "feature_id"]
        assert eligibility.loc[2, "display_identifier"] == "P3;P4"
        assert page.preparation["tables"]["imputed_cells"].empty
        assert page.parent_choice.currentData() == prep_a
        page._start_statistics(); wait_for_run(page)
        assert page.statistics is not None, page.status.text()
        stat_a = page.statistics["metadata"]["run_id"]
        baseline = page.statistics["tables"]["all_results"].copy()
        assert len(baseline) == len(page.preparation["tables"]["prepared_matrix"])
        assert baseline.loc[0, "gene_symbol"] == "GENE00"
        assert baseline.loc[0, "uniprot_accession"] == "P00000"
        page.search.setText("gene00")
        assert page.results_table.rowCount() == 1
        page.search.clear()
        assert len(set(baseline.df_residual)) > 1
        assert "ROW000005" not in baseline.feature_id.tolist()
        assert "ROW000005" in page.statistics["tables"]["qualitative_candidates"].feature_id.tolist()
        assert page.results_table.rowCount() > 0
        assert page.model_table.rowCount() > 0
        assert page.graph_choice.count() >= 8
        exports = Path(folder) / "exports"
        exports.mkdir()
        with patch.object(page, "_save_path", return_value=str(exports / "results.csv")):
            page.result_tabs.setCurrentIndex(1); page._export_table()
        assert (exports / "results.csv").is_file()
        with patch.object(page, "_save_path", return_value=str(exports / "results.xlsx")):
            page._export_workbook()
        assert (exports / "results.xlsx").is_file()
        with patch.object(page, "_save_path", return_value=str(exports / "volcano.pdf")):
            page.graph_choice.setCurrentIndex(0); page._export_graph()
        assert (exports / "volcano.pdf").is_file()
        page.results_table.selectRow(0)
        with patch.object(page, "_save_path", return_value=str(exports / "feature.csv")):
            page._export_feature()
        assert (exports / "feature.csv").is_file()
        page.imputation.setCurrentIndex(page.imputation.findData("MinProb"))
        page._start_preparation(); wait_for_run(page)
        assert page.preparation is not None, page.status.text()
        prep_b = page.preparation["metadata"]["run_id"]
        assert prep_b != prep_a
        assert len(page.preparation["tables"]["imputed_cells"]) > 0
        assert "ROW000005" not in page.preparation["tables"]["prepared_matrix"].feature_id.tolist()
        page.trend.setChecked(True); page.effect.setValue(1)
        page._start_statistics(); wait_for_run(page)
        assert page.statistics is not None, page.status.text()
        stat_b = page.statistics["metadata"]["run_id"]
        assert stat_b != stat_a
        assert page.statistics["metadata"]["parent_preparation_run_id"] == prep_b
        stat_b_results = page.statistics["tables"]["all_results"].copy()
        page.effect.setValue(100)
        page._start_statistics(); wait_for_run(page)
        stat_c_results = page.statistics["tables"]["all_results"]
        assert stat_c_results["adj.P.Val"].equals(stat_b_results["adj.P.Val"])
        assert (stat_c_results["combined_significant"] == "TRUE").sum() <= (stat_b_results["combined_significant"] == "TRUE").sum()
        page.imputation.setCurrentIndex(page.imputation.findData("none"))
        page.transformation.setProperty("user_changed", True)
        page.transformation.setCurrentIndex(page.transformation.findData("none"))
        page._start_preparation(); wait_for_run(page)
        assert page.scale.currentData() is None
        assert not page.statistics_button.isEnabled()
        page.scale.setCurrentIndex(page.scale.findData("continuous"))
        page._start_statistics(); wait_for_run(page)
        continuous = page.statistics["tables"]["all_results"]
        assert page.statistics["metadata"]["prepared_scale"] == "continuous"
        assert (continuous["log2FC"] == "").all()
        assert (continuous["fold_change"] == "").all()
        page.prep_history.setCurrentIndex(page.prep_history.findData(prep_a))
        assert page.preparation["metadata"]["run_id"] == prep_a
        page.prep_history.setCurrentIndex(page.prep_history.findData(prep_b))
        page.prep_history.setCurrentIndex(page.prep_history.findData(prep_a))
        page.stat_history.setCurrentIndex(page.stat_history.findData(stat_a))
        assert page.statistics["metadata"]["parent_preparation_run_id"] == prep_a
        page.stat_history.setCurrentIndex(page.stat_history.findData(stat_b))
        page.stat_history.setCurrentIndex(page.stat_history.findData(stat_a))
        assert page.statistics["metadata"]["parent_preparation_run_id"] == prep_a
        source.write_text("ProteinGroup,A1,A2,A3,B1,B2,B3\nChanged,1,1,1,1,1,1\n", encoding="utf-8")
        page._load_statistics_history()
        assert page.statistics["tables"]["all_results"].equals(baseline)
        parent_root = Path(project.root) / "analyses/Differential/Preparation/runs" / prep_a
        hidden = parent_root.with_name(parent_root.name + "_hidden")
        parent_root.rename(hidden)
        try:
            page._load_statistics_history()
            assert page.statistics["metadata"]["parent_preparation_run_id"] == prep_a
            assert page.parent_audit_tables["feature_eligibility"].rowCount() > 0
        finally:
            hidden.rename(parent_root)
        graph = page.statistics["run_root"] / "plots/volcano.png"
        graph.rename(graph.with_suffix(".hidden"))
        try:
            page._load_statistics_history()
            page.graph_choice.setCurrentIndex(0)
            assert page.graph_label.text() == "Graph not available for this run."
        finally:
            graph.with_suffix(".hidden").rename(graph)
        incomplete = Path(project.root) / "analyses/Differential/Statistics/runs/incomplete_run"
        incomplete.mkdir()
        page._refresh_statistics_history()
        page.stat_history.setCurrentIndex(page.stat_history.findData("incomplete_run"))
        assert page.statistics is None
        assert page.results_table.rowCount() == 0
        assert page.graph_choice.count() == 0
        print("Differential GUI offline smoke passed")
    page.close()


if __name__ == "__main__":
    smoke()
