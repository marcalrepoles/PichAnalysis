import hashlib
import json

import pytest

from pichanalysis.core.differential_analysis import (
    DifferentialParameters, NoReadyPreparationError, PreparationRunIncompleteError,
    PreparedMatrixHashMismatchError, PreparedScaleRequiredError,
    UnsupportedPreparedScaleError, MissingDifferentialOutputError,
    list_preparation_runs, list_runs, load_run, prepare_run, run_differential_analysis)
from pichanalysis.core.differential_preparation import (
    PreparationParameters, run_differential_preparation)
from pichanalysis.core.r_runtime import RRuntime
from tests.smoke_differential_preparation import fixture_project
from tests.smoke_differential_analysis import smoke


def parent(tmp_path, *, transformation=None):
    project, columns, source = fixture_project(tmp_path)
    output = run_differential_preparation(project, RRuntime(), run_id="parent",
        parameters=PreparationParameters("A", "B", tuple(columns), transformation=transformation))
    return project, output, source


def test_no_ready_parent(tmp_path):
    project, _, _ = fixture_project(tmp_path)
    with pytest.raises(NoReadyPreparationError):
        prepare_run(project, "statistics", DifferentialParameters("missing"))


def test_frozen_inputs_and_hashes(tmp_path):
    project, output, _ = parent(tmp_path)
    assert list_preparation_runs(project) == ["parent"]
    run, provenance = prepare_run(project, "statistics", DifferentialParameters("parent"))
    metadata = json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))
    assert metadata["prepared_scale"] == "log2"
    assert metadata["parent_preparation_run_id"] == "parent"
    assert metadata["prepared_matrix_hash"] == hashlib.sha256((run / "input/prepared_matrix.csv").read_bytes()).hexdigest()
    assert (provenance / "13_differential_analysis.R").is_file()
    assert (run / "input/qualitative_detection_candidates.csv").is_file()


def test_prepared_scale_required_for_untransformed_parent(tmp_path):
    project, _, _ = parent(tmp_path, transformation="none")
    with pytest.raises(PreparedScaleRequiredError):
        prepare_run(project, "bad", DifferentialParameters("parent"))
    run, _ = prepare_run(project, "continuous", DifferentialParameters("parent", prepared_scale="continuous"))
    assert json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))["effect_scale"] == "continuous"
    with pytest.raises(UnsupportedPreparedScaleError):
        prepare_run(project, "bad", DifferentialParameters("parent", prepared_scale="raw"))


def test_parent_hash_protection(tmp_path):
    project, output, _ = parent(tmp_path)
    matrix = output["run_root"] / "matrices/prepared_matrix.csv"
    matrix.write_text(matrix.read_text(encoding="utf-8").replace("ROW000001", "ROWXXXXXX"), encoding="utf-8")
    with pytest.raises(PreparedMatrixHashMismatchError):
        prepare_run(project, "statistics", DifferentialParameters("parent"))


def test_incomplete_parent(tmp_path):
    project, output, _ = parent(tmp_path)
    (output["run_root"] / "summary.csv").unlink()
    assert list_preparation_runs(project) == []
    with pytest.raises(PreparationRunIncompleteError):
        prepare_run(project, "statistics", DifferentialParameters("parent"))


def test_effect_threshold_does_not_change_bh(tmp_path):
    project, _, _ = parent(tmp_path)
    a = run_differential_analysis(project, RRuntime(), run_id="a",
        parameters=DifferentialParameters("parent", minimum_absolute_effect=0))
    b = run_differential_analysis(project, RRuntime(), run_id="b",
        parameters=DifferentialParameters("parent", minimum_absolute_effect=100))
    assert a["tables"]["all_results"]["adj.P.Val"].equals(b["tables"]["all_results"]["adj.P.Val"])
    assert len(b["tables"]["significant_results"]) == 0
    assert b["tables"]["significant_results"].columns.tolist() == b["tables"]["all_results"].columns.tolist()
    assert set(list_runs(project)) == {"a", "b"}
    assert load_run(project, "a")["metadata"]["minimum_absolute_effect"] == 0


def test_historical_missing_artifact(tmp_path):
    project, _, _ = parent(tmp_path)
    result = run_differential_analysis(project, RRuntime(), run_id="statistics",
        parameters=DifferentialParameters("parent"))
    (result["run_root"] / "results/differential_results_all.csv").unlink()
    with pytest.raises(MissingDifferentialOutputError):
        load_run(project, "statistics")


def test_offline_scientific_smoke():
    smoke()


def test_continuous_scale_end_to_end(tmp_path):
    project, _, _ = parent(tmp_path, transformation="none")
    output = run_differential_analysis(project, RRuntime(), run_id="continuous_stats",
        parameters=DifferentialParameters("parent", prepared_scale="continuous"))
    results = output["tables"]["all_results"]
    assert output["metadata"]["prepared_scale"] == "continuous"
    assert (results.effect_scale == "continuous").all()
    assert (results.log2FC == "").all()
    assert (results.fold_change == "").all()
    assert (results.percent_change == "").all()
    assert (results.difference_A_minus_B != "").any()