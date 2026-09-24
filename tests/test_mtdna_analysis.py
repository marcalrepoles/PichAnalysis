"""Offline gene-level mtDNA Evidence backend and history regressions."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from pichanalysis.core.mtdna_analysis import (
    MtdnaParameters, MtdnaTargetOutsideBackgroundError, MtdnaUnavailableError,
    MtdnaMissingOutputError, list_mtdna_runs, map_mtdna_experiment,
    prepare_mtdna_run, read_mtdna_outputs, run_mtdna_analysis,
)
from pichanalysis.core.r_runtime import RRuntime
from smoke_mtdna_analysis import build_fixture
from test_mtdna_evidence_database import sources


def test_readiness_and_gene_level_mapping(tmp_path):
    project, manager = build_fixture(tmp_path)
    mapping = map_mtdna_experiment(project, manager.mtdna_evidence)
    assert mapping.loc[mapping.source_row.eq("1"), "entity_key"].iloc[0] == "NCBI:1"
    assert mapping.loc[mapping.source_row.eq("21"), "entity_key"].iloc[0] == "NCBI:1"
    assert mapping.loc[mapping.source_row.eq("4"), "mapping_status"].iloc[0] == "mapped_no_mtdna_evidence"
    assert mapping.loc[mapping.source_row.eq("22"), "mapping_status"].iloc[0] == "ambiguous"
    assert mapping.loc[mapping.source_row.eq("23"), "entity_key"].iloc[0] == "NCBI:100"
    catalog_path = project.root / "mapping/tables/protein_catalog.csv"
    catalog = pd.read_csv(catalog_path, dtype=str, keep_default_na=False)
    catalog.loc[catalog.source_row.eq("1"), "uniprot_accession"] = "P22222"
    catalog.to_csv(catalog_path, index=False)
    conflicted = map_mtdna_experiment(project, manager.mtdna_evidence)
    assert conflicted.loc[conflicted.source_row.eq("1"), "mapping_status"].iloc[0] == "ambiguous"
    compatible = catalog.iloc[[0]].copy()
    compatible["source_row"] = "24"
    compatible["mapping_status"] = "ambiguous"
    compatible["uniprot_accession"] = "P11111-2;P11111-3"
    catalog = pd.concat([catalog, compatible], ignore_index=True)
    catalog.to_csv(catalog_path, index=False)
    resolved_group = map_mtdna_experiment(project, manager.mtdna_evidence)
    assert resolved_group.loc[resolved_group.source_row.eq("24"), "entity_key"].iloc[0] == "NCBI:1"
    manager.mtdna_evidence.active_pointer.unlink()
    with pytest.raises(MtdnaUnavailableError):prepare_mtdna_run(project,manager,"unavailable")


def test_experimental_background_dedup_safeguard_and_provenance(tmp_path):
    project, manager = build_fixture(tmp_path)
    with pytest.raises(MtdnaTargetOutsideBackgroundError) as error:
        prepare_mtdna_run(project, manager, "outside", MtdnaParameters(
            target_selection="Manual selection",manual_rows=(1,4),
            background_selection="Manual selection",background_manual_rows=(1,)))
    assert error.value.details == {"entities_outside_background":["NCBI:4"],
        "initial_target_size":2,"initial_background_size":1}
    run, provenance = prepare_mtdna_run(project, manager, "authorized", MtdnaParameters(
        target_selection="Manual selection",manual_rows=(1,4),
        background_selection="Manual selection",background_manual_rows=(1,),
        allow_target_outside_background=True))
    metadata = json.loads((run / "metadata.json").read_text())
    assert metadata["target_size"] == metadata["background_size"] == 1
    assert metadata["target_adjustment_authorized"]
    assert (provenance / "mtdna_evidence_manifest.json").is_file()
    assert (provenance / "10_mtdna_analysis.R").is_file()
    run2, _ = prepare_mtdna_run(project, manager, "all")
    assert len(pd.read_csv(run2 / "inputs/background.csv")) == 21
    assert len(pd.read_csv(run2 / "inputs/target.csv")) == 21
    assert "NCBI:4" in set(pd.read_csv(run2 / "inputs/background.csv").entity_key)
    no_evidence = run_mtdna_analysis(project, manager, RRuntime(), run_id="no-evidence",
        parameters=MtdnaParameters(target_selection="Manual selection", manual_rows=(4,),
            minimum_overlap=1))
    assert no_evidence["metadata"]["background_size"] == 21
    assert no_evidence["tables"]["source_count"].loc[0,"entity_count"] == "1"
    assert no_evidence["tables"]["enrichment_significant"].empty


def test_r_outputs_history_snapshot_lock_and_missing_artifact(tmp_path):
    project, manager = build_fixture(tmp_path)
    runtime = RRuntime()
    a = run_mtdna_analysis(project,manager,runtime,run_id="run-A",
        parameters=MtdnaParameters(target_selection="A-specific",minimum_overlap=2))
    b = run_mtdna_analysis(project,manager,runtime,run_id="run-B",
        parameters=MtdnaParameters(target_selection="B-specific",minimum_overlap=1))
    assert a["metadata"]["target_size"] == 3
    assert a["metadata"]["background_size"] == 21
    assert len(a["tables"]["enrichment_significant"]) >= 1
    assert a["workbook"].is_file() and len(a["plots"]) == 6
    assert b["metadata"]["target_definition"] == "B-specific"
    assert read_mtdna_outputs(project,"run-A")["metadata"]["run_id"] == "run-A"
    assert list_mtdna_runs(project) == ["run-B","run-A"]
    snapshot_x = a["metadata"]["snapshot_id"]
    mito, go = sources(manager.root,"fixture-B")
    manager.mtdna_evidence.build(mitocarta=mito,go_snapshot=go,ncbi_xml=tmp_path / "ncbi.xml")
    y = run_mtdna_analysis(project,manager,runtime,run_id="run-Y",
        parameters=MtdnaParameters(target_selection="A-specific",minimum_overlap=1))
    assert y["metadata"]["snapshot_id"] != snapshot_x
    assert read_mtdna_outputs(project,"run-A")["metadata"]["snapshot_id"] == snapshot_x
    (y["run_root"] / "plots/category_frequency.png").unlink()
    assert len(read_mtdna_outputs(project,"run-Y")["plots"]) == 5
    (y["run_root"] / "summary.csv").unlink()
    with pytest.raises(MtdnaMissingOutputError):read_mtdna_outputs(project,"run-Y")
