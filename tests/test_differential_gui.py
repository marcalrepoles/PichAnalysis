"""Fast, offline checks for Differential Analysis GUI decisions."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from pichanalysis.ui.differential_analysis_page import DifferentialAnalysisPage
from tests.smoke_differential_preparation import fixture_project


@pytest.fixture
def page(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, _, _ = fixture_project(tmp_path)
    widget = DifferentialAnalysisPage()
    widget.set_project(project)
    yield widget
    widget.close()


def test_condition_and_sample_validation(page):
    assert page.prepare_button.isEnabled()
    assert page.selected_columns() == ("A1", "A2", "A3", "B1", "B2", "B3")
    page.condition_b.setCurrentIndex(page.condition_a.currentIndex())
    assert not page.prepare_button.isEnabled()
    page.condition_b.setCurrentIndex(1)
    page.samples.item(0, 0).setCheckState(Qt.CheckState.Unchecked)
    page.samples.item(1, 0).setCheckState(Qt.CheckState.Unchecked)
    assert not page.prepare_button.isEnabled()


def test_defaults_and_explicit_options(page):
    assert page.transformation.currentData() == "log2_positive"
    assert page.zero_missing.isChecked()
    assert page.normalization.currentData() == "none"
    assert page.minimum_observed.value() == 2
    assert page.imputation.currentData() == "none"
    assert page.seed.value() == 12345
    params = page.preparation_parameters()
    assert params.imputation_method == "none"
    assert params.condition_a == "A" and params.condition_b == "B"
    page.normalization.setCurrentIndex(page.normalization.findData("median_center"))
    page.imputation.setCurrentIndex(page.imputation.findData("MinProb"))
    page.q.setValue(.02); page.sigma.setValue(1.2); page.margin.setCurrentIndex(0)
    params = page.preparation_parameters()
    assert (params.normalization, params.imputation_method, params.minprob_q,
        params.imputation_sigma, params.imputation_margin) == ("median_center", "MinProb", .02, 1.2, 1)
    page.imputation.setCurrentIndex(page.imputation.findData("KNN"))
    page.k.setValue(4)
    assert page.preparation_parameters().knn_k == 4


def test_spectral_counts_blocked(page):
    for row in range(page.samples.rowCount()):
        page.samples.item(row, 0).setData(Qt.ItemDataRole.UserRole, "spectral_count")
    page._update_state()
    assert not page.prepare_button.isEnabled()
    assert "count-based model" in page.preparation_readiness.text()


def test_other_continuous_requires_confirmation(page):
    for row in range(page.samples.rowCount()):
        page.samples.item(row, 0).setData(Qt.ItemDataRole.UserRole, "other_quantitative")
    page._update_state()
    assert not page.prepare_button.isEnabled()
    page.other_confirmed.setChecked(True)
    assert page.prepare_button.isEnabled()


def test_statistics_defaults(page):
    assert page.scale.currentData() is None
    assert not page.trend.isChecked()
    assert not page.robust.isChecked()
    assert page.fdr.value() == .05
    assert page.effect.value() == 0
    assert page.top_n.value() == 20
    assert page.volcano_labels.value() == 20


def test_no_worker_without_ready_parent(page):
    page._start_statistics()
    assert not page.is_running()
