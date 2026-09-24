import pytest

from pichanalysis.core.differential_preparation import (
    PreparationParameters, prepare_run, run_differential_preparation, load_run, list_runs,
    TooFewConditionsError, SameConditionError, TooFewSamplesError,
    MixedQuantificationTypesError, SpectralCountsUnsupportedError,
    OtherContinuousUnconfirmedError, MissingSampleMetadataError,
    InvalidNumericValueError, NegativeValueError, ZeroLog2IncompatibleError,
    UnsupportedTransformationError, UnsupportedNormalizationError,
    UnsupportedImputationError, NoEligibleFeaturesError)
from pichanalysis.core.r_runtime import RRuntime
from tests.smoke_differential_preparation import fixture_project


def parameters(columns, **kwargs):
    return PreparationParameters("A", "B", tuple(columns), **kwargs)


def test_freezes_rows_and_metadata(tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    run, provenance = prepare_run(project, "freeze", parameters(columns))
    assert (run / "input/original_quantitative_matrix.csv").is_file()
    assert (provenance / "12_differential_preparation.R").is_file()
    assert len((run / "input/feature_metadata.csv").read_text(encoding="utf-8").splitlines()) == 13


def test_conditions_and_samples(tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    with pytest.raises(SameConditionError):
        prepare_run(project, "bad", PreparationParameters("A", "A", tuple(columns)))
    with pytest.raises(TooFewSamplesError):
        prepare_run(project, "bad", parameters(("A1", "B1", "B2")))
    with pytest.raises(TooFewConditionsError):
        prepare_run(project, "bad", PreparationParameters("A", "C", tuple(columns)))
    project.config["columns"]["A1"]["replicate"] = ""
    with pytest.raises(MissingSampleMetadataError):
        prepare_run(project, "bad", parameters(columns))


def test_family_validation(tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    project.config["columns"]["B1"]["quantification_type"] = "raw_intensity"
    with pytest.raises(MixedQuantificationTypesError):
        prepare_run(project, "bad", parameters(columns))
    for col in columns:
        project.config["columns"][col]["quantification_type"] = "spectral_count"
    with pytest.raises(SpectralCountsUnsupportedError):
        prepare_run(project, "bad", parameters(columns))
    for col in columns:
        project.config["columns"][col]["quantification_type"] = "other_quantitative"
    with pytest.raises(OtherContinuousUnconfirmedError):
        prepare_run(project, "bad", parameters(columns))
    prepare_run(project, "other", parameters(columns, other_continuous_confirmed=True,
        zero_is_missing=True, transformation="none"))


@pytest.mark.parametrize("field,value,error", [
    ("transformation", "log2p1", UnsupportedTransformationError),
    ("normalization", "quantile", UnsupportedNormalizationError),
    ("imputation_method", "minimum", UnsupportedImputationError),
])
def test_unsupported_options(tmp_path, field, value, error):
    project, columns, _ = fixture_project(tmp_path)
    with pytest.raises(error):
        prepare_run(project, "bad", parameters(columns, **{field: value}))


def test_invalid_values_and_log2_safeguards(tmp_path):
    project, columns, path = fixture_project(tmp_path)
    original = path.read_text(encoding="utf-8")
    path.write_text(original.replace("P1,100,110", "P1,typo,110"), encoding="utf-8")
    with pytest.raises(InvalidNumericValueError) as error:
        prepare_run(project, "bad", parameters(columns))
    assert error.value.details["invalid_value_count"] == 1
    path.write_text(original.replace("P1,100,110", "P1,-100,110"), encoding="utf-8")
    with pytest.raises(NegativeValueError):
        prepare_run(project, "bad", parameters(columns))
    prepare_run(project, "negative_none", parameters(columns, transformation="none"))
    path.write_text(original, encoding="utf-8")
    with pytest.raises(ZeroLog2IncompatibleError):
        prepare_run(project, "bad", parameters(columns, zero_is_missing=False))
    prepare_run(project, "zero_none", parameters(columns, zero_is_missing=False, transformation="none"))


def test_no_eligible_features_and_qualitative_audit(tmp_path):
    project, columns, path = fixture_project(tmp_path)
    path.write_text("ProteinGroup," + ",".join(columns) + "\nAonly,10,11,12,,,\nBonly,,,,20,21,22\n", encoding="utf-8")
    with pytest.raises(NoEligibleFeaturesError):
        run_differential_preparation(project, RRuntime(), run_id="qualitative_only",
            parameters=parameters(columns))
    run = project.root / "analyses/Differential/Preparation/runs/qualitative_only"
    assert (run / "eligibility/qualitative_detection_candidates.csv").is_file()
    assert "qualitative_only" in list_runs(project)
    from pichanalysis.core.differential_preparation import MissingPreparationOutputError
    with pytest.raises(MissingPreparationOutputError):
        load_run(project, "qualitative_only")


@pytest.mark.parametrize("method", ["MinProb", "QRILC", "KNN"])
def test_official_imputation_methods(tmp_path, method):
    project, columns, _ = fixture_project(tmp_path)
    result = run_differential_preparation(project, RRuntime(), run_id=method,
        parameters=parameters(columns, imputation_method=method, knn_k=3))
    assert result["metadata"]["imputation_method"] == method
    assert len(result["tables"]["imputed_cells"]) > 0
    assert "ROW000005" not in result["tables"]["prepared_matrix"].feature_id.tolist()


def test_structured_r_dependency_and_execution_errors(tmp_path):
    from types import SimpleNamespace
    from pichanalysis.core.differential_preparation import ImputationPackageUnavailableError, PreparationRExecutionError
    project, columns, _ = fixture_project(tmp_path)
    class Runtime:
        def __init__(self, stderr): self.stderr = stderr
        def run(self, *args, **kwargs): return SimpleNamespace(returncode=1, stderr=self.stderr, stdout="")
    with pytest.raises(ImputationPackageUnavailableError):
        run_differential_preparation(project, Runtime("IMPUTATION_PACKAGE_UNAVAILABLE: MsCoreUtils"),
            run_id="dependency_failure", parameters=parameters(columns, imputation_method="MinProb"))
    with pytest.raises(PreparationRExecutionError):
        run_differential_preparation(project, Runtime("synthetic R failure"),
            run_id="r_failure", parameters=parameters(columns))

def test_nonfinite_and_zero_counts_are_frozen(tmp_path):
    import json
    project, columns, path = fixture_project(tmp_path)
    path.write_text(path.read_text(encoding="utf-8").replace("P1,100,110", "P1,Inf,NA"), encoding="utf-8")
    run, _ = prepare_run(project, "counts", parameters(columns))
    metadata = json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))
    assert metadata["input_value_counts"]["Inf"] == 1
    assert metadata["input_value_counts"]["NA"] == 1
    assert metadata["input_value_counts"]["zero"] == 1

def test_run_id_stays_within_project(tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    from pichanalysis.core.differential_preparation import DifferentialPreparationError
    with pytest.raises(DifferentialPreparationError):
        prepare_run(project, "../outside", parameters(columns))