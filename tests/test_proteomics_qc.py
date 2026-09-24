import pytest

from pichanalysis.core.proteomics_qc import (
    QCParameters, NoQuantitativeColumnsError, TooFewSamplesError,
    MixedQuantificationTypesError, MissingConditionError, MissingReplicateError,
    InvalidNumericValueError, NegativeValueError, UnsupportedTransformationError,
    NoUsableQuantitativeDataError, prepare_qc_run,
)
from tests.smoke_proteomics_qc import fixture_project


def test_prepare_freezes_original_rows_and_hashes(tmp_path):
    project, cols, _ = fixture_project(tmp_path)
    run, provenance = prepare_qc_run(project, "unit", QCParameters(tuple(cols)))
    assert (run / "input/feature_metadata.csv").is_file()
    assert (provenance / "proteomics_qc.R").is_file()


@pytest.mark.parametrize("selected,error", [
    ((), NoQuantitativeColumnsError),
    (("A1",), TooFewSamplesError),
    (("A1", "bad"), NoQuantitativeColumnsError),
])
def test_sample_selection_errors(tmp_path, selected, error):
    project, _, _ = fixture_project(tmp_path)
    with pytest.raises(error):
        prepare_qc_run(project, "invalid", QCParameters(selected))


def test_mixed_family_rejected(tmp_path):
    project, cols, _ = fixture_project(tmp_path)
    project.config["columns"]["B1"]["quantification_type"] = "spectral_count"
    with pytest.raises(MixedQuantificationTypesError):
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols)))


@pytest.mark.parametrize("field,error", [("condition", MissingConditionError), ("replicate", MissingReplicateError)])
def test_required_metadata(tmp_path, field, error):
    project, cols, _ = fixture_project(tmp_path)
    project.config["columns"]["A1"][field] = ""
    with pytest.raises(error):
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols)))


def test_invalid_numeric_value(tmp_path):
    project, cols, path = fixture_project(tmp_path)
    contents = path.read_text(encoding="utf-8").replace("P1,10,12", "P1,typo,12")
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(InvalidNumericValueError) as caught:
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols)))
    assert caught.value.details["invalid_value_count"] == 1


def test_negative_intensity_rejected(tmp_path):
    project, cols, path = fixture_project(tmp_path)
    path.write_text(path.read_text(encoding="utf-8").replace("P1,10,12", "P1,-10,12"), encoding="utf-8")
    with pytest.raises(NegativeValueError):
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols)))


def test_unsupported_transform_rejected(tmp_path):
    project, cols, _ = fixture_project(tmp_path)
    with pytest.raises(UnsupportedTransformationError):
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols), transformation="sqrt"))


def test_all_zero_intensity_rejected(tmp_path):
    project, cols, path = fixture_project(tmp_path)
    path.write_text("ProteinGroup," + ",".join(cols) + "\nP," + ",".join(["0"] * len(cols)) + "\n", encoding="utf-8")
    with pytest.raises(NoUsableQuantitativeDataError):
        prepare_qc_run(project, "invalid", QCParameters(tuple(cols)))
