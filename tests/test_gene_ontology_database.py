"""Offline contract tests for managed human GO snapshots."""
from __future__ import annotations

import gzip
import json

import pytest

from pichanalysis.core.gene_ontology_database import (
    GeneOntologyCancelled, GeneOntologyDatabase, GeneOntologyError, parse_gaf,
)


def _sources(path):
    path.mkdir()
    (path / "go-basic.obo").write_text(
        "format-version: 1.2\ndata-version: releases/2026-01-01\n\n"
        "[Term]\nid: GO:0000001\nname: parent\nnamespace: biological_process\n\n"
        "[Term]\nid: GO:0000002\nname: child\nnamespace: biological_process\n"
        "is_a: GO:0000001 ! parent\n\n"
        "[Term]\nid: GO:0000003\nname: part\nnamespace: biological_process\n"
        "relationship: part_of GO:0000002 ! child\n", encoding="utf-8")
    fields = ["UniProtKB", "P12345", "TFAM", "", "GO:0000003", "PMID:1", "IDA", "", "P",
              "", "", "protein", "taxon:9606", "20260101", "UniProt", "", ""]
    excluded = fields.copy()
    excluded[3] = "NOT"
    with gzip.open(path / "HUMAN-uniprot.gaf.gz", "wt", encoding="utf-8") as stream:
        stream.write("!gaf-version: 2.2\n")
        stream.write("\t".join(fields) + "\n")
        stream.write("\t".join(excluded) + "\n")


def test_go_snapshot_queries_hashes_and_update_rollback(tmp_path):
    sources = tmp_path / "sources"
    _sources(sources)
    db = GeneOntologyDatabase(tmp_path / "databases")
    first = db.build_from_directory(sources)
    assert db.active_snapshot() == first
    assert db.validate_snapshot(first)
    assert db.manifest(first)["annotation_version"] == "2.2"
    assert db.manifest(first)["counts"]["not_annotations"] == 1
    assert db.get_descendants("GO:0000001").go_id.tolist() == [
        "GO:0000001", "GO:0000002", "GO:0000003"]
    annotations = db.get_annotations_for_terms(["GO:0000003"])
    assert len(annotations) == 2
    assert set(annotations.qualifier) == {"", "NOT"}
    assert set(annotations.evidence_code) == {"IDA"}
    assert set(annotations.reference) == {"PMID:1"}
    assert set(annotations.assigned_by) == {"UniProt"}
    with pytest.raises(GeneOntologyCancelled):
        db.build_from_directory(sources, cancel_requested=lambda: True)
    assert db.active_snapshot() == first
    (sources / "go-basic.obo").write_text("bad", encoding="utf-8")
    with pytest.raises(GeneOntologyError):
        db.build_from_directory(sources)
    assert db.active_snapshot() == first
    _sources_repair(sources)
    assert db.build_from_directory(sources) != first
    assert db.active_snapshot() != first
    assert db.validate_snapshot(first)
    (first / "raw" / "go-basic.obo").write_text("corrupted", encoding="utf-8")
    with pytest.raises(GeneOntologyError):
        db.validate_snapshot(first)


def _sources_repair(path):
    path.joinpath("go-basic.obo").write_text(
        "[Term]\nid: GO:0000001\nname: parent\n"
        "[Term]\nid: GO:0000002\nname: child\nis_a: GO:0000001 ! parent\n",
        encoding="utf-8")

def test_gaf_rejects_nonhuman(tmp_path):
    sources = tmp_path / "sources"
    _sources(sources)
    gaf = sources / "HUMAN-uniprot.gaf.gz"
    with gzip.open(gaf, "wt", encoding="utf-8") as stream:
        stream.write("\t".join(["UniProtKB", "P12345", "TFAM", "", "GO:0000003",
            "PMID:1", "IDA", "", "P", "", "", "protein", "taxon:10090", "20260101",
            "UniProt", "", ""]) + "\n")
    with pytest.raises(GeneOntologyError, match="no annotations"):
        parse_gaf(gaf)
