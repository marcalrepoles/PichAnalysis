"""Offline Differential -> QC / Presence Explorer smoke with frozen inputs."""
import os
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from pichanalysis.core.cross_module_integration import (open_cross_module_index,
    rebuild_cross_module_index)
from pichanalysis.core.differential_analysis import DifferentialParameters, run_differential_analysis
from pichanalysis.core.differential_preparation import PreparationParameters, run_differential_preparation
from pichanalysis.core.presence_analysis import build_presence_arguments
from pichanalysis.core.proteomics_qc import QCParameters, run_proteomics_qc
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.analyses_page import AnalysesPage
from smoke_differential_preparation import fixture_project
from test_cross_module_adapters import test_multi_module_offline_bridges_from_persisted_ids
from test_cross_module_navigation import test_legacy_mapping_never_borrows_latest


def smoke():
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory() as folder:
        project, columns, source = fixture_project(Path(folder))
        runtime = RRuntime()
        run_proteomics_qc(project, runtime, run_id="qc_a",
            parameters=QCParameters(tuple(columns)))
        args = build_presence_arguments(project, quantification_type="lfq_intensity",
            selected_conditions=["A", "B"], zero_is_missing=True, threshold=0,
            rule_mode="count", rule_value=2, predominant=False, run_id="presence_a")
        result = runtime.run(Path(__file__).resolve().parents[1] / "r_scripts/02_presence_absence.R",
            *args, timeout=180)
        assert result.returncode == 0, result.stderr
        run_differential_preparation(project, runtime, run_id="prep_a",
            parameters=PreparationParameters("A", "B", tuple(columns)))
        run_differential_analysis(project, runtime, run_id="stat_a",
            parameters=DifferentialParameters("prep_a"))
        index = rebuild_cross_module_index(project)
        assert index.healthy()
        sources = [row for row in index.search("ROW000001")
            if row["module"] == "differential"]
        assert len(sources) == 1
        context = index.context(sources[0]["id"])
        assert context.feature_id == "ROW000001" and context.source_row == 1
        assert context.uniprot_accessions == ()
        matches = index.matches(context)
        assert any(match.target_module == "proteomics_qc" and match.status.value == "Exact"
            and match.lineage_compatible for match in matches)
        assert any(match.target_module == "presence_absence" and match.status.value == "Exact"
            and match.lineage_compatible for match in matches), [(m.target_module, m.status.value, m.basis) for m in matches]
        page = AnalysesPage()
        page.set_project(project)
        page.show()
        page.explore_context("differential", "stat_a", "ROW000001")
        assert page.module_tabs.currentWidget() is page.cross_module_explorer
        assert page.cross_module_explorer.context.feature_id == "ROW000001"
        assert page.cross_module_explorer.matches_table.rowCount() >= 2
        page.open_analysis_target("proteomics_qc", "qc_a", "ROW000001")
        assert page.module_tabs.currentWidget() is page.proteomics_qc_page
        assert page.proteomics_qc_page.outputs["metadata"]["run_id"] == "qc_a"
        assert "Target: ROW000001" in page.proteomics_qc_page.cross_module_target_hint.text()
        assert page.proteomics_qc_page.tables[("Detection", "feature_detection")].currentRow() >= 0, page.proteomics_qc_page.cross_module_target_hint.text()
        page.proteomics_qc_page.explore_button.click()
        assert page.cross_module_explorer.context.source_module == "proteomics_qc"
        assert page.cross_module_explorer.context.source_run_id == "qc_a"
        page.open_analysis_target("presence_absence", "presence_a", "P1", source_row="1")
        assert page.module_tabs.currentWidget() is page.presence_page
        assert page.presence_page.outputs.metadata["run_id"] == "presence_a"
        assert "Target: P1" in page.presence_page.cross_module_target_hint.text()
        assert page.presence_page.table.currentRow() >= 0
        page.presence_page.cross_module_explore_button.click()
        assert page.cross_module_explorer.context is not None, (page.presence_page.cross_module_explore_hint.text(), page.cross_module_explorer.status.text(), page.cross_module_explorer.search.text(), page.cross_module_explorer.source_records)
        assert page.cross_module_explorer.context.source_module == "presence_absence"
        assert page.cross_module_explorer.context.source_run_id == "presence_a"
        assert any(match.target_module == "proteomics_qc" and match.status.value == "Exact"
            for match in page.cross_module_explorer.matches)
        page.open_analysis_target("differential", "stat_a", "ROW000001")
        assert page.module_tabs.currentWidget() is page.differential_analysis_page
        assert page.differential_analysis_page.statistics["metadata"]["run_id"] == "stat_a"
        page.differential_analysis_page.results_table.selectRow(0)
        page.differential_analysis_page.explore_button.click()
        assert page.cross_module_explorer.context.source_module == "differential"
        assert page.cross_module_explorer.context.source_run_id == "stat_a"
        source.write_text("ProteinGroup,A1,A2,A3,B1,B2,B3\nChanged,1,1,1,1,1,1\n", encoding="utf-8")
        again = open_cross_module_index(project)
        assert again.context(sources[0]["id"]).original_identifier == "P1"
        index.path.write_bytes(b"corrupt index")
        recovered = open_cross_module_index(project)
        assert recovered.healthy()
        assert len(recovered.search("ROW000001")) >= 2
        page.close()
        print("Cross-module integration offline GUI smoke passed")


if __name__ == "__main__":
    smoke()
