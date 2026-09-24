"""Offline Python -> R -> Qt GUI smoke for Proteomics QC."""
import os
import tempfile
from pathlib import Path

import pandas as pd

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QTimer, Qt
from PySide6.QtWidgets import QApplication, QFileDialog

from pichanalysis.ui.main_window import MainWindow
from tests.smoke_proteomics_qc import fixture_project


def wait_for_worker(page, timeout_ms=30000):
    worker = page._workers[-1]
    loop = QEventLoop()
    timed_out = []
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(lambda: (timed_out.append(True), loop.quit()))
    worker.finished.connect(loop.quit)
    timer.start(timeout_ms)
    loop.exec()
    timer.stop()
    assert not timed_out, "Proteomics QC worker timed out"
    QApplication.processEvents()
    assert not page.is_running(), page.status.text()
    assert page.outputs is not None, page.status.text()


def smoke():
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory() as folder:
        project, columns, source = fixture_project(Path(folder))
        data = pd.read_csv(source)
        data.loc[data.ProteinGroup == "P6", columns] = 0
        data.loc[data.ProteinGroup == "P2", "B3"] = float("nan")
        data.to_csv(source, index=False)
        window = MainWindow()
        page = window.analyses_page.proteomics_qc_page
        window.analyses_page.set_project(project)
        window.navigation.setCurrentRow(2)
        window.show()
        assert page.samples.rowCount() == 6
        assert page.selected_columns() == tuple(columns)
        assert page.transformation.currentData() == "log2_positive"
        assert page.zero_missing.isChecked()
        assert page.run_button.isEnabled()
        page._start_run()
        assert not page.run_button.isEnabled()
        wait_for_worker(page)
        a = page.outputs
        a_id = a["metadata"]["run_id"]
        a_summary = a["tables"]["summary"].copy()
        a_graphs = [page.graph_choice.itemText(i) for i in range(page.graph_choice.count())]
        assert len(a_graphs) == 11
        assert page.graph_preview.pixmap() is not None
        assert len(page._frames["sample_metrics"]) == 6
        for key in ("sample_detection", "feature_detection", "sample_missingness",
            "pairwise_correlations", "detection_overlap", "feature_condition_cv",
            "pca_scores", "sample_diagnostics", "warnings"):
            assert key in page._frames
        assert len(page._frames["feature_detection"]) == 8
        assert page._frames["feature_detection"].loc[0, "feature_id"] != page._frames["feature_detection"].loc[1, "feature_id"]
        assert page._frames["pca_summary"].loc[0, "status"] == "Ready"
        assert int(page._frames["pca_summary"].loc[0, "complete_case_features_used"]) == 2
        assert page._frames["transformed_matrix"].loc[6, "S001"] == ""
        assert (page._frames["pairwise_correlations"].pearson == "").any()
        assert any(table.item(row, col).text() == "NA" for table in [page.tables[("Correlations", "pairwise_correlations")]] for row in range(table.rowCount()) for col in range(table.columnCount()))
        assert "quality_score" not in page._frames["sample_diagnostics"]
        assert abs(float(page._frames["feature_condition_cv"].query("feature_id == 'ROW000001' and condition == 'A'").iloc[0].cv_percent) - 100*2/12) < 1e-8
        assert "score" not in page.history_metadata.text().lower()
        page.transformation.setCurrentIndex(page.transformation.findData("none"))
        page.zero_missing.setChecked(False)
        page._start_run()
        wait_for_worker(page)
        b_id = page.outputs["metadata"]["run_id"]
        assert b_id != a_id and page.outputs["metadata"]["transformation"] == "none"
        source.write_text("ProteinGroup,A1,A2,A3,B1,B2,B3\nChanged,1,1,1,1,1,1\n", encoding="utf-8")
        page.history.setCurrentIndex(page.history.findData(a_id))
        assert page.outputs["metadata"]["run_id"] == a_id
        assert page.outputs["metadata"]["transformation"] == "log2_positive"
        assert page.outputs["metadata"]["zero_is_missing"] is True
        assert page._frames["summary"].equals(a_summary)
        assert len(page._frames["sample_metadata"]) == 6
        assert [page.graph_choice.itemText(i) for i in range(page.graph_choice.count())] == a_graphs
        old_frame = page._frames["feature_detection"].copy()
        page.tabs.setCurrentIndex(list(page.table_tabs).index("Detection"))
        nested = page.table_tabs["Detection"]
        nested.setCurrentIndex(1)
        page.searches[("Detection", "feature_detection")].setText("P1")
        assert page.tables[("Detection", "feature_detection")].rowCount() == 2
        assert page._frames["feature_detection"].equals(old_frame)
        destination = Path(folder) / "export.csv"
        original_dialog = QFileDialog.getSaveFileName
        QFileDialog.getSaveFileName = staticmethod(lambda *args, **kwargs: (str(destination), ""))
        try:
            page._export_current()
            assert len(destination.read_text(encoding="utf-8").splitlines()) == 3
            destination = Path(folder) / "export.xlsx"
            page._export_workbook()
            assert destination.is_file()
            destination = Path(folder) / "export.png"
            page._export_graph()
            assert destination.is_file()
            destination = Path(folder) / "export.pdf"
            page._export_graph()
            assert destination.is_file()
        finally:
            QFileDialog.getSaveFileName = original_dialog
        window.close()
        assert not window.isVisible()
        print("Proteomics QC Python -> R -> GUI smoke passed")


if __name__ == "__main__":
    smoke()
