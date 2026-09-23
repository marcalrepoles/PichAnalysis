from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.complex_analysis import (
    ComplexParameters, ComplexPortalUnavailableError, ComplexPresenceRequiredError,
    ComplexRExecutionError, ComplexTargetOutsideBackgroundError,
    ComplexUnsupportedOrganismError, InvalidComplexComponentError,
    MissingComplexOutputError, NoComplexProteinsError, list_complex_runs,
    map_complex_experiment, prepare_complex_run, read_complex_outputs,
    run_complex_analysis,
)
from pichanalysis.core.r_runtime import RResult, RRuntime
from smoke_complex_analysis import build_fixture


def test_mapping_canonical_target_background_and_safeguard(tmp_path):
    project, manager = build_fixture(tmp_path)
    mapping = map_complex_experiment(project, manager.complex_portal)
    assert mapping.loc[mapping.source_row == "1", "canonical_uniprot"].iloc[0] == "P12345"
    assert mapping.loc[mapping.source_row == "2", "canonical_uniprot"].iloc[0] == "P12345"
    assert mapping.loc[mapping.source_row == "2", "isoform_normalized"].iloc[0]
    assert mapping.loc[mapping.source_row == "10", "membership_status"].iloc[0] == "no_curated_complex_membership"
    assert mapping.loc[mapping.source_row == "11", "mapping_status"].iloc[0] == "ambiguous"
    assert mapping.loc[mapping.source_row == "12", "mapping_status"].iloc[0] == "unmapped"
    with pytest.raises(ComplexTargetOutsideBackgroundError) as error:
        prepare_complex_run(project, manager, "outside", ComplexParameters(
            target_selection="Manual selection", manual_rows=(1,3),
            background_selection="Manual selection", background_manual_rows=(1,)))
    assert error.value.details == {"entities_outside_background": ["Q11111"],
        "initial_target_size": 2, "initial_background_size": 1}
    run, provenance = prepare_complex_run(project, manager, "authorized", ComplexParameters(
        target_selection="Manual selection", manual_rows=(1,3),
        background_selection="Manual selection", background_manual_rows=(1,),
        allow_target_outside_background=True))
    assert pd.read_csv(run / "inputs/target.csv").canonical_uniprot.tolist() == ["P12345"]
    assert str(pd.read_csv(run / "inputs/target.csv").target_supporting_rows.iloc[0]) == "1"
    assert (provenance / "complex_portal_manifest.json").is_file()


def test_presence_and_structured_prerequisite_errors(tmp_path):
    project, manager = build_fixture(tmp_path)
    with pytest.raises(ComplexPresenceRequiredError):
        prepare_complex_run(project, manager, "presence", ComplexParameters(target_selection="Shared"))
    classes = project.root / "analyses/presence_absence/tables/classification.csv"
    classes.parent.mkdir(parents=True)
    pd.DataFrame([{"source_row": "1", "classification": "Shared"},
                  {"source_row": "3", "classification": "Condition-specific"}]).to_csv(classes, index=False)
    run, _ = prepare_complex_run(project, manager, "shared", ComplexParameters(target_selection="Shared"))
    assert pd.read_csv(run / "inputs/target.csv").canonical_uniprot.tolist() == ["P12345"]
    assert str(pd.read_csv(run / "inputs/target.csv").target_supporting_rows.iloc[0]) == "1"
    project.config["organism_tax_id"] = "10090"
    with pytest.raises(ComplexUnsupportedOrganismError):
        prepare_complex_run(project, manager, "mouse")
    project.config["organism_tax_id"] = "9606"
    manager.complex_portal.active_pointer.unlink()
    with pytest.raises(ComplexPortalUnavailableError):
        prepare_complex_run(project, manager, "missing_database")


def test_no_unique_proteins_and_r_error_mapping(tmp_path):
    project, manager = build_fixture(tmp_path)
    class FailedRuntime:
        def run(self, *args, **kwargs):
            return RResult((), 1, "", "Invalid Complex Portal component representation: invalid")
    with pytest.raises(InvalidComplexComponentError):
        run_complex_analysis(project, manager, FailedRuntime(), run_id="invalid_component")
    class CrashedRuntime:
        def run(self, *args, **kwargs):
            raise RuntimeError("R unavailable")
    with pytest.raises(ComplexRExecutionError):
        run_complex_analysis(project, manager, CrashedRuntime(), run_id="r_failure")
    catalog = project.root / "mapping/tables/protein_catalog.csv"
    pd.DataFrame([{"source_row": 1, "original_id": "x", "uniprot_accession": "NONE"}]).to_csv(catalog, index=False)
    with pytest.raises(NoComplexProteinsError):
        prepare_complex_run(project, manager, "no_proteins")


def test_r_outputs_history_and_missing_artifact(tmp_path):
    project, manager = build_fixture(tmp_path)
    runtime = RRuntime()
    if not runtime.available:
        pytest.skip("Rscript unavailable")
    params = ComplexParameters(target_selection="Manual selection", manual_rows=(1,2,3,4,5,6))
    a = run_complex_analysis(project, manager, runtime, run_id="a", parameters=params, timeout=180)
    first_snapshot = a["metadata"]["snapshot_id"]
    assert a["tables"]["summary"].shape[0] >= 20
    assert len(a["plots"]) == 12
    assert a["workbook"].is_file()
    second = manager.complex_portal.install_from_file(tmp_path / "9606.tsv")
    b = run_complex_analysis(project, manager, runtime, run_id="b", parameters=ComplexParameters(
        target_selection="Manual selection", manual_rows=(1,3)), timeout=180)
    assert second.name == b["metadata"]["snapshot_id"] != first_snapshot
    historical = read_complex_outputs(project, "a")
    assert historical["metadata"]["snapshot_id"] == first_snapshot
    assert historical["tables"]["complex_coverage"].equals(a["tables"]["complex_coverage"])
    assert list_complex_runs(project) == ["a", "b"]
    (a["run_root"] / "coverage/component_groups.csv").unlink()
    with pytest.raises(MissingComplexOutputError):
        read_complex_outputs(project, "a")
