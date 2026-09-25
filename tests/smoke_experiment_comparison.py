"""Offline GUI smoke for two tables, two frozen Differential runs, and history."""
import hashlib
import json
import os
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
from PySide6.QtWidgets import QApplication

from pichanalysis.core.differential_analysis import TABLES
from pichanalysis.ui.experiment_comparison_page import ExperimentComparisonPage


def frozen_differential(project_root: Path, run_id: str, directions: list[str]):
    root = project_root / "analyses/Differential/Statistics/runs" / run_id
    (root / "input").mkdir(parents=True)
    prepared = root / "input/prepared_matrix.csv"
    prepared.write_text("feature_id,value\nF1,1\n", encoding="utf-8")
    metadata = {"run_id": run_id, "prepared_matrix_hash":
        hashlib.sha256(prepared.read_bytes()).hexdigest(),
        "condition_A": "Induced", "condition_B": "Control",
        "comparison_direction": "Condition A - Condition B", "effect_scale": "log2"}
    (root / "input/parameters.json").write_text(json.dumps(metadata), encoding="utf-8")
    for name in ("sample_metadata.csv", "feature_metadata.csv", "contrast_metadata.csv",
                 "preparation_summary.csv", "preparation_parameters.json", "feature_eligibility.csv",
                 "qualitative_detection_candidates.csv"):
        (root / "input" / name).write_text("item,value\nplaceholder,1\n", encoding="utf-8")
    results = pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12347"],
        "effect_direction": directions, "combined_significant": ["TRUE", "TRUE", "FALSE"],
        "effect": ["1", "-1", "0.2"], "log2FC": ["1", "-1", "0.2"],
        "P.Value": ["0.001", "0.001", "0.4"], "adj.P.Val": ["0.01", "0.01", "0.5"]})
    for key, relative in TABLES.items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        (results if key == "all_results" else pd.DataFrame({"placeholder": []})).to_csv(destination, index=False)
    (root / "Differential_analysis.xlsx").write_bytes(b"frozen-smoke-fixture")


def smoke():
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory() as folder:
        root = Path(folder)
        a, b = root / "dataset_a.csv", root / "dataset_b.csv"
        pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12347", "A99999"],
            "description": ["one", "two", "three", "four"]}).to_csv(a, index=False)
        pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12347", "B99999"],
            "description": ["one", "two", "three", "five"]}).to_csv(b, index=False)
        frozen_differential(root, "diff_a", ["higher_in_condition_A", "higher_in_condition_B",
            "higher_in_condition_A"])
        frozen_differential(root, "diff_b", ["higher_in_condition_A", "higher_in_condition_A",
            "higher_in_condition_B"])
        project = SimpleNamespace(root=root, config={"input": {"processed_file": None}})
        page = ExperimentComparisonPage()
        page.set_project(project)
        page.source_a.set_source(a)
        page.source_b.set_source(b)
        for source in (page.source_a, page.source_b):
            source.column.setCurrentText("uniprot_accession")
            source.kind.setCurrentIndex(source.kind.findData("uniprot"))
            source.confirmed.setChecked(True)
        page.comparison_type.setCurrentIndex(page.comparison_type.findData("uniprot"))
        page.differential_a.setCurrentIndex(page.differential_a.findData("diff_a"))
        page.differential_b.setCurrentIndex(page.differential_b.findData("diff_b"))
        assert "Induced" in page.contrasts.text()
        page.match_a.setChecked(True)
        page.match_b.setChecked(True)
        page.comparable.setChecked(True)
        page._compare()
        assert page.loaded is not None, page.status.text()
        assert len(page.loaded["tables"]["a_only"]) == 1
        assert len(page.loaded["tables"]["b_only"]) == 1
        assert len(page.loaded["tables"]["shared"]) == 3
        assert "Significant opposite direction" in set(page.loaded["tables"]
            ["differential_comparison"].comparison_class)
        export = root / "export.xlsx"
        with patch("pichanalysis.ui.experiment_comparison_page.QFileDialog.getSaveFileName",
                   return_value=(str(export), "")):
            page._export_workbook()
        with pd.ExcelFile(export) as workbook:
            assert {"Summary", "Master", "Differential", "A only", "B only", "Shared"} <= set(
                workbook.sheet_names)
        run_id = page.loaded["run_root"].name
        a.unlink()
        b.unlink()
        page.set_project(project)
        page.history.setCurrentIndex(page.history.findData(run_id))
        page._load_history()
        assert page.loaded["run_root"].name == run_id
        page.close()
    print("Experiment Comparison offline GUI smoke: passed")


if __name__ == "__main__":
    smoke()
