import json
from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.r_runtime import RResult
from pichanalysis.core.reactome_analysis import (
    MissingReactomeOutputError, ReactomeAnalysisError, ReactomeDatabaseUnavailableError, ReactomeRExecutionError,
    ReactomeParameters, TargetOutsideBackgroundError, UnsupportedReactomeOrganismError,
    list_reactome_runs, prepare_reactome_arguments, reactome_readiness,
    read_reactome_outputs, run_reactome_analysis, select_reactome_set,
)


def reactome_source(root: Path) -> Path:
    root.mkdir()
    files = {
        "ReactomePathways.txt": (
            "R-HSA-1\tPathway one\tHomo sapiens\n"
            "R-HSA-2\tPathway two\tHomo sapiens\n"
            "R-HSA-3\tShared child\tHomo sapiens\n"
            "R-MMU-1\tMouse\tMus musculus\n"
        ),
        "ReactomePathwaysRelation.txt": "R-HSA-1\tR-HSA-3\nR-HSA-2\tR-HSA-3\nR-MMU-1\tR-MMU-2\n",
        "UniProt2Reactome.txt": (
            "U1\tR-HSA-1\thttps://x\tPathway one\tIEA\tHomo sapiens\n"
            "U1\tR-HSA-2\thttps://x\tPathway two\tIEA\tHomo sapiens\n"
            "U1\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens\n"
            "U2\tR-HSA-1\thttps://x\tPathway one\tIEA\tHomo sapiens\n"
            "U2\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens\n"
            "M1\tR-MMU-1\thttps://x\tMouse\tIEA\tMus musculus\n"
        ),
        "NCBI2Reactome.txt": (
            "N1\tR-HSA-2\thttps://x\tPathway two\tIEA\tHomo sapiens\n"
            "N1\tR-HSA-3\thttps://x\tShared child\tIEA\tHomo sapiens\n"
        ),
        "humanPathwaysWithDiagrams.txt": "R-HSA-1\nR-HSA-2\nR-HSA-3\n",
        "pathway2summation.txt": "R-HSA-1\tSummary one\nR-HSA-2\tSummary two\nR-HSA-3\tSummary three\n",
    }
    for name, contents in files.items():
        (root / name).write_text(contents, encoding="utf-8")
    return root


def reactome_project(tmp_path, *, organism="Homo sapiens", tax_id="9606"):
    project = create_project(tmp_path, "ReactomeProject")
    set_organism(project, organism, tax_id)
    catalog = pd.DataFrame([
        {"source_row": 1, "original_id": "row1", "mapping_status": "mapped_unique", "uniprot_accession": "U1", "ncbi_gene_id": "N1", "gene_symbol": "G1", "protein_name": "Protein 1"},
        {"source_row": 2, "original_id": "row2", "mapping_status": "mapped_unique", "uniprot_accession": None, "ncbi_gene_id": "N1", "gene_symbol": "G2", "protein_name": "Protein 2"},
        {"source_row": 3, "original_id": "row3", "mapping_status": "mapped_unique", "uniprot_accession": "U2", "ncbi_gene_id": None, "gene_symbol": "G3", "protein_name": "Protein 3"},
        {"source_row": 4, "original_id": "U1;U2", "mapping_status": "ambiguous", "uniprot_accession": "U1;U2", "ncbi_gene_id": None, "gene_symbol": "GROUP", "protein_name": "Protein group"},
        {"source_row": 5, "original_id": "unknown", "mapping_status": "unmapped", "uniprot_accession": None, "ncbi_gene_id": None, "gene_symbol": "X", "protein_name": "Unknown"},
    ])
    catalog.to_csv(project.root / "mapping" / "tables" / "protein_catalog.csv", index=False)
    classification = pd.DataFrame({
        "source_row": [1, 2, 3, 4, 5],
        "classification": ["A-specific", "Sporadic", "Shared", "Shared", "Not reproducibly detected"],
    })
    tables = project.root / "analyses" / "presence_absence" / "tables"
    tables.mkdir(parents=True)
    classification.to_csv(tables / "classification.csv", index=False)
    return project


def manager_with_snapshot(tmp_path):
    manager = DatabaseManager(tmp_path / "databases")
    snapshot = manager.reactome.install_from_directory(reactome_source(tmp_path / "reactome-source"))
    manifest = manager.reactome.manifest(snapshot)
    manifest["release_version"] = "97"
    manager.reactome._write_manifest(snapshot, manifest)
    return manager, snapshot


def test_parameters_pin_active_snapshot_and_core_paths(tmp_path):
    project = reactome_project(tmp_path)
    manager, snapshot = manager_with_snapshot(tmp_path)
    args = prepare_reactome_arguments(project, manager, run_id="parameters", parameters=ReactomeParameters(minimum_overlap=2))
    assert args[args.index("--snapshot-id") + 1] == snapshot.name
    assert args[args.index("--minimum-overlap") + 1] == "2"
    hashes = json.loads(args[args.index("--snapshot-hashes") + 1])
    assert len(hashes) == 6
    assert Path(args[args.index("--uniprot-map") + 1]).parent == snapshot / "raw"


def test_missing_database_is_structured(tmp_path):
    project = reactome_project(tmp_path)
    manager = DatabaseManager(tmp_path / "empty-db")
    assert not reactome_readiness(project, manager).ready
    with pytest.raises(ReactomeDatabaseUnavailableError, match="Reactome database unavailable"):
        prepare_reactome_arguments(project, manager, run_id="missing")


def test_unsupported_organism_is_structured(tmp_path):
    project = reactome_project(tmp_path, organism="Mus musculus", tax_id="10090")
    manager, _ = manager_with_snapshot(tmp_path)
    with pytest.raises(UnsupportedReactomeOrganismError, match="Unsupported organism"):
        prepare_reactome_arguments(project, manager, run_id="mouse")


def test_presence_and_manual_target_selection(tmp_path):
    project = reactome_project(tmp_path)
    assert set(select_reactome_set(project, "Shared").source_row) == {3, 4}
    assert list(select_reactome_set(project, "Manual selection", (2,)).source_row) == [2]
    (project.root / "analyses" / "presence_absence" / "tables" / "classification.csv").unlink()
    with pytest.raises(ReactomeAnalysisError, match="Presence/Absence outputs"):
        select_reactome_set(project, "Shared")


def test_missing_expected_output_is_structured(tmp_path):
    project = reactome_project(tmp_path)
    with pytest.raises(MissingReactomeOutputError, match="Missing expected output"):
        read_reactome_outputs(project, "absent")


def test_r_execution_failure_is_structured(tmp_path):
    project = reactome_project(tmp_path)
    manager, _ = manager_with_snapshot(tmp_path)

    class FailedRuntime:
        def run(self, script, *arguments, timeout):
            return RResult((str(script),), 1, "", "planned R failure")

    with pytest.raises(ReactomeRExecutionError, match="R execution failure: planned R failure"):
        run_reactome_analysis(project, manager, FailedRuntime(), run_id="r-failure")


def test_target_outside_background_and_authorized_adjustment(tmp_path):
    project = reactome_project(tmp_path)
    manager, _ = manager_with_snapshot(tmp_path)
    runtime = RRuntime()
    assert runtime.available
    blocked = ReactomeParameters(target_selection="Shared", background_selection="A-specific", minimum_overlap=1)
    with pytest.raises(TargetOutsideBackgroundError) as caught:
        run_reactome_analysis(project, manager, runtime, run_id="blocked", parameters=blocked)
    assert caught.value.details["initial_target_size"] == 1
    assert caught.value.details["initial_background_size"] == 1
    assert caught.value.details["entities_outside_background"] == ["UP:U2"]

    allowed = ReactomeParameters(target_selection="Shared", background_selection="A-specific",
        minimum_overlap=1, allow_target_outside_background=True)
    outputs = run_reactome_analysis(project, manager, runtime, run_id="allowed", parameters=allowed)
    adjustment = outputs.metadata["target_background_adjustment"]
    assert adjustment["user_decision"] == "Continue"
    assert adjustment["final_target_size"] == 0


def test_python_to_r_scientific_smoke_and_persistence(tmp_path):
    project = reactome_project(tmp_path)
    manager, snapshot = manager_with_snapshot(tmp_path)
    runtime = RRuntime()
    assert runtime.available
    outputs = run_reactome_analysis(project, manager, runtime, run_id="smoke",
        parameters=ReactomeParameters(minimum_overlap=1, fdr_cutoff=0.05, top_n=20))
    status = outputs.mapping.drop_duplicates("source_row")["mapping_status"].value_counts().to_dict()
    assert status == {"mapped_unique": 3, "ambiguous": 1, "unmapped": 1}
    assert len(outputs.membership) == 7
    assert len(outputs.hierarchy) == 2
    assert len(outputs.frequency) == 3
    assert len(outputs.enrichment) == 3
    assert outputs.workbook.is_file() and outputs.workbook.stat().st_size > 0
    assert len(outputs.graphs) == 6 and all(path.stat().st_size > 0 for path in outputs.graphs)
    assert outputs.metadata["snapshot_id"] == snapshot.name
    assert outputs.metadata["source"] == "Local Reactome snapshot"
    assert outputs.run_root.is_dir() and outputs.provenance_root.is_dir()
    assert (outputs.run_root / "raw" / "target_entities.csv").is_file()
    assert (outputs.provenance_root / "parameters.json").is_file()
    assert (outputs.provenance_root / "database_hashes.json").is_file()
    assert (outputs.provenance_root / "05_reactome_analysis.R").is_file()
    assert list_reactome_runs(project) == ["smoke"]
    latest = read_reactome_outputs(project)
    assert latest.metadata["run_id"] == "smoke"
    summary_values = dict(zip(outputs.summary["metric"], outputs.summary["value"].astype(str)))
    print("REACTOME_SCIENTIFIC_SMOKE " + json.dumps({
        "input_entities": 5,
        "mapped_unique": status.get("mapped_unique", 0),
        "ambiguous": status.get("ambiguous", 0),
        "unmapped": status.get("unmapped", 0),
        "pathway_memberships": len(outputs.membership),
        "hierarchy_rows": len(outputs.hierarchy),
        "frequency_rows": len(outputs.frequency),
        "eligible_pathways": int(float(summary_values["Pathways eligible"])),
        "tested_pathways": len(outputs.enrichment),
        "significant_pathways": len(outputs.significant),
        "plots": len(outputs.graphs),
        "workbook": str(outputs.workbook),
        "run_directory": str(outputs.run_root),
        "provenance_directory": str(outputs.provenance_root),
    }, sort_keys=True))
