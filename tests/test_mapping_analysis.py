import json
from pathlib import Path

import pandas as pd

import pichanalysis.core.mapping_analysis as module
from pichanalysis.core.mapping_analysis import (
    build_mapping_arguments, cache_path, export_result, mapping_readiness, read_mapping_outputs,
)
from pichanalysis.core.organism import get_organism, set_organism
from pichanalysis.core.project import create_project, open_project


def configured_project(tmp_path):
    project = create_project(tmp_path, "Analysis")
    processed = project.root / "input" / "processed" / "data.csv"
    processed.write_text("Protein IDs\nP04637\n", encoding="utf-8")
    project.config["input"]["processed_file"] = "input/processed/data.csv"
    project.config["input"]["rows"] = 1
    project.config["columns"] = {"Protein IDs": {
        "role": "identifier", "identifier_type": "uniprot", "primary_identifier": True,
        "condition": "", "replicate": "", "quantification_type": "unknown", "detection": {},
    }}
    project.save()
    return project


def test_organism_persists(tmp_path):
    project = create_project(tmp_path, "Organism")
    set_organism(project, "Homo sapiens", "9606")
    reopened = open_project(project.root)
    assert get_organism(reopened).name == "Homo sapiens"
    assert get_organism(reopened).tax_id == "9606"


def test_analysis_blocked_without_primary_identifier(tmp_path):
    project = create_project(tmp_path, "NoIdentifier")
    set_organism(project, "Homo sapiens", "9606")
    state = mapping_readiness(project)
    assert not state.ready
    assert "identificador principal" in state.reason.lower()


def test_analysis_blocked_without_organism(tmp_path):
    state = mapping_readiness(configured_project(tmp_path))
    assert not state.ready
    assert "organismo" in state.reason.lower()


def test_safe_mapping_arguments(tmp_path):
    project = configured_project(tmp_path)
    set_organism(project, "Homo sapiens", "9606")
    arguments = build_mapping_arguments(project, tmp_path / "cache.sqlite", False, "run-1")
    assert arguments[arguments.index("--id-column") + 1] == "Protein IDs"
    assert arguments[arguments.index("--tax-id") + 1] == "9606"
    assert arguments[arguments.index("--refresh") + 1] == "false"
    assert all(isinstance(value, str) for value in arguments)


def test_cache_uses_platform_location(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "user_cache_path", lambda *args: tmp_path / "cache root")
    result = cache_path()
    assert result == tmp_path / "cache root" / "annotation_cache.sqlite"
    assert result.parent.is_dir()


def test_reads_r_results(tmp_path):
    project = create_project(tmp_path, "Results")
    tables = project.root / "mapping" / "tables"
    frame = pd.DataFrame({"input_id": ["P04637"], "mapping_status": ["mapped_unique"]})
    for name in ("protein_catalog.csv", "id_mapping.csv", "unmapped.csv", "ambiguous.csv"):
        frame.to_csv(tables / name, index=False)
    (project.root / "mapping" / "latest_metadata.json").write_text(
        json.dumps({"mapped_unique_count": 1}), encoding="utf-8")
    outputs = read_mapping_outputs(project)
    assert outputs.catalog.loc[0, "input_id"] == "P04637"
    assert outputs.metadata["mapped_unique_count"] == 1


def test_export_does_not_change_original(tmp_path):
    source = tmp_path / "source.csv"
    destination = tmp_path / "elsewhere" / "copy.csv"
    source.write_bytes(b"id\nP04637\n")
    before = source.read_bytes()
    export_result(source, destination)
    assert source.read_bytes() == before
    assert destination.read_bytes() == before
