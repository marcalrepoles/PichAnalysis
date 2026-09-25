import os
from types import SimpleNamespace

import pandas as pd
from PySide6.QtWidgets import QApplication

from pichanalysis.ui.experiment_comparison_page import ExperimentComparisonPage


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_identity_options_require_safe_resolution(tmp_path):
    app = QApplication.instance() or QApplication([])
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    pd.DataFrame({"id": ["P12345"], "description": ["A"]}).to_csv(a, index=False)
    pd.DataFrame({"id": ["GENE1"], "description": ["B"]}).to_csv(b, index=False)
    page = ExperimentComparisonPage()
    page.set_project(SimpleNamespace(root=tmp_path, config={"input": {"processed_file": None}}))
    page.source_a.set_source(a)
    page.source_b.set_source(b)
    page.source_a.kind.setCurrentIndex(page.source_a.kind.findData("uniprot"))
    page.source_b.kind.setCurrentIndex(page.source_b.kind.findData("gene_symbol"))
    assert page.comparison_type.count() == 0
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame({"uniprot_accession": ["P12345"], "gene_symbol": ["GENE1"]}).to_csv(
        mapping, index=False)
    page.mapping_path = mapping
    page._refresh_identity_options()
    assert {page.comparison_type.itemData(i) for i in range(page.comparison_type.count())} == {
        "uniprot", "gene_symbol"}
    page.close()
