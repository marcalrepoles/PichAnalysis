"""Offline end-to-end Proteomics QC smoke, including frozen historical runs."""
import hashlib
import tempfile
from pathlib import Path

import pandas as pd

from pichanalysis.core.project import create_project
from pichanalysis.core.proteomics_qc import QCParameters, TABLES, PLOTS, load_run, run_proteomics_qc
from pichanalysis.core.r_runtime import RRuntime


def fixture_project(root):
    project = create_project(root, "QCFixture")
    path = project.root / "input/processed/qc.csv"
    columns = ["A1", "A2", "A3", "B1", "B2", "B3"]
    rows = [
        ["P1", 10, 12, 14, 20, 22, 24],
        ["P1", 20, 24, 28, 40, 44, 48],
        ["P2", 40, 42, 44, 80, 82, 84],
        ["P3", 0, 5, 0, 7, 0, 8],
        ["P4", 10, 0, 0, 0, 0, 0],
        ["P5", 0, 0, 0, 50, 55, 60],
        ["P6", 5, 5, 5, 5, 5, 5],
        ["P7", "", "", "", "", "", ""],
    ]
    pd.DataFrame(rows, columns=["ProteinGroup", *columns]).to_csv(path, index=False)
    project.config["input"].update({"processed_file": "input/processed/qc.csv", "rows": len(rows)})
    project.config["columns"] = {"ProteinGroup": {"role": "identifier", "identifier_type": "protein_group", "primary_identifier": True}}
    for col in columns:
        project.config["columns"][col] = {"role": "quantification", "condition": col[0],
            "replicate": col[1], "quantification_type": "lfq_intensity"}
    project.save()
    return project, columns, path


def smoke():
    with tempfile.TemporaryDirectory() as directory:
        project, columns, source = fixture_project(Path(directory))
        runtime = RRuntime()
        a = run_proteomics_qc(project, runtime, run_id="smoke_a",
            parameters=QCParameters(tuple(columns)))
        assert set(TABLES) == set(a["tables"])
        assert a["workbook"].is_file()
        assert len(a["plots"]) == len(PLOTS)
        for image in a["plots"]:
            assert image.stat().st_size > 1000
            assert image.with_suffix(".pdf").is_file()
        assert len(a["tables"]["feature_detection"]) == 8
        assert a["tables"]["feature_detection"]["display_identifier"].tolist()[:2] == ["P1", "P1"]
        assert a["metadata"]["feature_unit"] == "one original experimental row"
        for name, digest in a["metadata"]["input_sha256"].items():
            assert hashlib.sha256((a["run_root"] / "input" / name).read_bytes()).hexdigest() == digest
        # Scientific invariants, checked against known fixture values.
        assert a["tables"]["transformed_matrix"].loc[3, "S001"] == ""
        assert a["tables"]["detection_matrix"].loc[3, "S001"] == "0"
        assert a["tables"]["feature_detection"].loc[0, "feature_id"] != a["tables"]["feature_detection"].loc[1, "feature_id"]
        cv = a["tables"]["feature_condition_cv"]
        known_cv = cv[(cv.feature_id == "ROW000001") & (cv.condition == "A")].iloc[0]
        assert abs(float(known_cv.cv_percent) - 100 * 2 / 12) < 1e-8
        pca = a["tables"]["pca_summary"].iloc[0]
        assert int(pca.complete_case_features_used) == 3
        assert pca.center == "TRUE" and pca.scale == "FALSE"
        pairs = a["tables"]["pairwise_correlations"]
        known_pair = pairs[(pairs.sample_a == "S001") & (pairs.sample_b == "S002")].iloc[0]
        assert int(known_pair.n_shared_features) == 4
        assert known_pair.pearson and known_pair.spearman
        before = (a["run_root"] / "input/quantitative_matrix.csv").read_bytes()
        b = run_proteomics_qc(project, runtime, run_id="smoke_b",
            parameters=QCParameters(tuple(columns), transformation="none", zero_is_missing=False))
        assert b["metadata"]["transformation"] == "none"
        source.write_text("ProteinGroup,A1,A2,A3,B1,B2,B3\nChanged,1,1,1,1,1,1\n", encoding="utf-8")
        old = load_run(project, "smoke_a")
        assert (old["run_root"] / "input/quantitative_matrix.csv").read_bytes() == before
        assert len(old["tables"]["feature_detection"]) == 8
        print("Proteomics QC offline smoke passed")


if __name__ == "__main__":
    smoke()
