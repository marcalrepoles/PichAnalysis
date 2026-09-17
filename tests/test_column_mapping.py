import pandas as pd

from pichanalysis.core.column_mapping import (
    ColumnRole,
    experimental_design,
    load_mapping,
    save_mapping,
    suggest_column,
    suggest_columns,
    validate_mapping,
)
from pichanalysis.core.project import create_project, open_project


def test_role_suggestions():
    identifier = suggest_column("Protein IDs", pd.Series(["P12345", "Q99999"]))
    quantitative = suggest_column("LFQ intensity Control 1", pd.Series([1.2, 3.4]))
    annotation = suggest_column("Description", pd.Series(["mitochondrial protein"]))
    unknown = suggest_column("Notes", pd.Series(["first sample", "second sample"]))
    assert identifier.role == ColumnRole.IDENTIFIER.value
    assert quantitative.role == ColumnRole.QUANTIFICATION.value
    assert quantitative.condition == "Control"
    assert quantitative.replicate == "1"
    assert annotation.role == ColumnRole.ANNOTATION.value
    assert unknown.role == ColumnRole.OTHER.value


def valid_columns():
    return {
        "Protein IDs": {
            "role": "identifier", "identifier_type": "uniprot", "primary_identifier": True,
            "condition": "", "replicate": "", "quantification_type": "unknown", "detection": {},
        },
        "LFQ Control 1": {
            "role": "quantification", "identifier_type": "unknown", "primary_identifier": False,
            "condition": "Control", "replicate": "1", "quantification_type": "lfq_intensity", "detection": {},
        },
    }


def test_valid_primary_identifier_and_design_api(tmp_path):
    project = create_project(tmp_path, "Mapping")
    columns = valid_columns()
    validation = save_mapping(project, columns)
    assert validation.valid
    design = experimental_design(project)
    assert design["primary_identifier_column"] == "Protein IDs"
    assert design["primary_identifier_type"] == "uniprot"
    assert design["conditions"] == ["Control"]


def test_mapping_persists_after_reopen(tmp_path):
    project = create_project(tmp_path, "Persistencia")
    columns = valid_columns()
    save_mapping(project, columns)
    reopened = open_project(project.root)
    assert load_mapping(reopened) == columns


def test_validation_without_identifier():
    result = validate_mapping({"Description": {"role": "annotation"}})
    assert not result.valid
    assert "No primary identifier selected." in result.errors


def test_validation_quantification_without_condition():
    columns = valid_columns()
    columns["LFQ Control 1"]["condition"] = ""
    result = validate_mapping(columns)
    assert any("has no condition" in error for error in result.errors)


def test_validation_duplicate_condition_replicate():
    columns = valid_columns()
    columns["LFQ Control copy"] = dict(columns["LFQ Control 1"])
    result = validate_mapping(columns)
    assert any("mesma combinação" in error for error in result.errors)


def test_suggestions_choose_one_primary_identifier():
    frame = pd.DataFrame({"Protein IDs": ["P12345"], "Gene names": ["TFAM"]})
    suggestions = suggest_columns(frame)
    assert sum(item["primary_identifier"] for item in suggestions.values()) == 1
