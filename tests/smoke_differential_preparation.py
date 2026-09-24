"""Offline Python -> R differential matrix-preparation smoke, without model fitting."""
import tempfile
from pathlib import Path

import pandas as pd

from pichanalysis.core.differential_preparation import (
    PreparationParameters, load_run, run_differential_preparation, TABLES, PLOTS)
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime


def fixture_project(root):
    project = create_project(root, "DifferentialFixture")
    columns = ["A1", "A2", "A3", "B1", "B2", "B3"]
    rows = [
        ["P1", 100, 110, 120, 200, 210, 220],
        ["P1", 150, 160, 170, 300, 310, 320],
        ["P3;P4", 80, 90, 100, 160, 170, 180],
        ["P5", 10, 11, "", 8, 9, ""],
        ["P6", 20, 21, 22, "", "", ""],
        ["P7", "", "", "", 30, 31, 32],
        ["P8", 12, "", "", "", "", ""],
        ["P9", "", "", "", "", "", ""],
        ["P10", 0, 15, 16, 25, 26, 27],
        ["P11", 40, 41, 42, 60, 61, 62],
        ["P12", 55, 56, 57, 65, 66, 67],
        ["P13", 70, 71, 72, 90, 91, 92],
    ]
    path = project.root / "input/processed/differential.csv"
    pd.DataFrame(rows, columns=["ProteinGroup", *columns]).to_csv(path, index=False)
    project.config["input"].update({"processed_file": "input/processed/differential.csv", "rows": len(rows)})
    project.config["columns"] = {"ProteinGroup": {"role": "identifier", "identifier_type": "protein_group", "primary_identifier": True}}
    for col in columns:
        project.config["columns"][col] = {"role": "quantification", "condition": col[0],
            "replicate": col[1], "quantification_type": "lfq_intensity"}
    project.save()
    return project, columns, path


def smoke():
    with tempfile.TemporaryDirectory() as folder:
        project, columns, path = fixture_project(Path(folder))
        runtime = RRuntime()
        a = run_differential_preparation(project, runtime, run_id="smoke_a",
            parameters=PreparationParameters("A", "B", tuple(columns)))
        assert set(a["tables"]) == set(TABLES)
        assert len(a["plots"]) == len(PLOTS)
        assert a["workbook"].is_file()
        eligibility = a["tables"]["feature_eligibility"]
        assert eligibility.loc[0, "pattern"] == "complete_both_conditions"
        assert eligibility.loc[3, "pattern"] == "partial_both_conditions"
        assert eligibility.loc[4, "pattern"] == "detected_only_condition_A"
        assert eligibility.loc[5, "pattern"] == "detected_only_condition_B"
        assert eligibility.loc[6, "pattern"] == "sparse_only_condition_A"
        assert eligibility.loc[7, "pattern"] == "all_missing"
        assert eligibility.loc[0, "feature_id"] != eligibility.loc[1, "feature_id"]
        assert eligibility.loc[2, "display_identifier"] == "P3;P4"
        prepared_a = a["tables"]["prepared_matrix"]
        assert len(prepared_a) < len(eligibility)
        assert prepared_a.loc[prepared_a.feature_id == "ROW000004", "S003"].iloc[0] == ""
        assert "ROW000005" not in prepared_a.feature_id.tolist()
        assert len(a["tables"]["qualitative_detection_candidates"]) == 2
        b = run_differential_preparation(project, runtime, run_id="smoke_b",
            parameters=PreparationParameters("A", "B", tuple(columns), normalization="median_center",
                imputation_method="MinProb", imputation_seed=12345))
        assert b["metadata"]["imputation_seed"] == 12345
        assert len(b["tables"]["imputed_cells"]) > 0
        assert "ROW000005" not in b["tables"]["prepared_matrix"].feature_id.tolist()
        normalized = b["tables"]["normalized_matrix"].set_index("feature_id")
        prepared = b["tables"]["prepared_matrix"].set_index("feature_id")
        for feature_id in prepared.index:
            for sample in [f"S{i:03d}" for i in range(1, 7)]:
                before = normalized.loc[feature_id, sample]
                after = prepared.loc[feature_id, sample]
                if before != "":
                    assert float(before) == float(after)
        path.write_text("ProteinGroup,A1,A2,A3,B1,B2,B3\nChanged,1,1,1,1,1,1\n", encoding="utf-8")
        old = load_run(project, "smoke_a")
        assert old["tables"]["prepared_matrix"].equals(prepared_a)
        assert old["metadata"]["normalization"] == "none"
        assert old["metadata"]["imputation_method"] == "none"
        print("Differential preparation offline smoke passed")


if __name__ == "__main__":
    smoke()
