import csv
import json
import sqlite3

import pytest

from pichanalysis.core.complex_portal_database import (
    ComplexPortalDatabase, canonicalize_uniprot_for_complex_portal,
)
from pichanalysis.core.databases.complex_portal import ComplexPortalCancelled

HEADERS = [
    "#Complex ac", "Recommended name", "Aliases for complex", "Taxonomy identifier",
    "Identifiers (and stoichiometry) of molecules in complex", "Evidence Code",
    "Experimental evidence", "Go Annotations", "Cross references", "Description",
    "Complex properties", "Complex assembly", "Ligand", "Disease", "Agonist",
    "Antagonist", "Comment", "Source", "Expanded participant list",
]


def fixture(path, *, tax="9606", predicted=False, conflict=False):
    rows = [
        ["CPX-1.2", "A complex", "Alpha|First", tax,
         "P12345(1)|CPX-2(0)|CHEBI:123(2)|URS000001(0)",
         "ECO:0008004(machine learning method evidence used in automatic assertion)" if predicted else "ECO:0000353(physical interaction evidence used in manual assertion)",
         "intact:EBI-1", "-", "-", "Description", "Property", "Assembly", "-", "-", "-", "-", "-",
         'psi-mi:"MI:0469"(IntAct)',
         "P12345(1)|Q11111(2)|Q22222(4)" + ("|Q11111(3)" if conflict else "")],
        ["CPX-2", "Nested complex", "Second", tax,
         "Q11111(2)|Q22222(4)", "ECO:0005547(manual assertion)", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",
         'psi-mi:"MI:0469"(IntAct)', "Q11111(2)|Q22222(4)"],
        ["CPX-3", "Another complex", "Third", tax,
         "P12345(0)", "ECO:0005547(manual assertion)", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",
         'psi-mi:"MI:0469"(IntAct)', "P12345(0)"],
    ]
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(HEADERS)
        writer.writerows(rows)
    return path


def test_curated_snapshot_queries_and_stoichiometry(tmp_path):
    source = fixture(tmp_path / "9606.tsv")
    database = ComplexPortalDatabase(tmp_path / "db")
    snapshot = database.install_from_file(source)
    manifest = database.manifest()
    assert database.active_snapshot() == snapshot
    assert database.validate_snapshot()
    assert manifest["scope"] == "manually_curated"
    assert manifest["complex_portal_release"] == "unknown"
    assert manifest["detected_headers"] == HEADERS
    assert manifest["complex_count"] == 3
    assert manifest["protein_membership_count"] == 6
    assert manifest["nonprotein_participant_count"] == 2
    assert manifest["nested_complex_count"] == 1
    assert database.source_hashes()["9606.tsv"] == manifest["raw_file"]["sha256"]
    item = database.get_complex("CPX-1")
    assert item["complex_accession_raw"] == "CPX-1.2"
    assert item["complex_version"] == "2"
    assert item["recommended_name"] == "A complex"
    assert item["source"] == 'psi-mi:"MI:0469"(IntAct)'
    assert item["confidence_accession"] == "ECO:0000353"
    assert item["confidence_name"] == "physical interaction evidence used in manual assertion"
    assert len(database.get_direct_participants("CPX-1")) == 4
    assert {row["participant_type"] for row in database.get_direct_participants("CPX-1")} == {"protein", "complex", "chemical", "RNA"}
    assert len(database.get_proteins_for_complex("CPX-1")) == 3
    assert database.get_nested_complexes("CPX-1")[0]["child_complex_id"] == "CPX-2"
    assert {row["participant_type"] for row in database.get_nonprotein_participants("CPX-1")} == {"chemical", "RNA"}
    unknown = database.get_nested_complexes("CPX-1")[0]
    assert unknown["stoichiometry_raw"] == "0" and unknown["stoichiometry"] == "" and unknown["stoichiometry_known"] == "0"
    assert {row["complex_id"] for row in database.get_complexes_for_uniprot("P12345")} == {"CPX-1", "CPX-3"}
    hits = database.get_complexes_for_uniprot("P12345-2")
    assert len(hits) == 2
    assert hits[0]["input_accession"] == "P12345-2"
    assert hits[0]["lookup_accession"] == "P12345"
    assert hits[0]["isoform_normalized"] is True
    assert database.get_complexes_for_uniprot("P12345")[0]["isoform_normalized"] is False
    assert database.search_complexes("alpha")[0]["complex_id"] == "CPX-1"
    assert database.search_complexes("q22222")
    assert database.search_complexes("cpx-1")
    database.provider = None
    assert database.validate_snapshot() and database.get_complex("CPX-2")


@pytest.mark.parametrize("value", ["P12345-bad", "P12345-0", "NOT-2", "CHEBI:123-2"])
def test_invalid_isoform_not_normalized(value):
    result = canonicalize_uniprot_for_complex_portal(value)
    assert result["lookup_accession"] == value
    assert not result["isoform_normalized"]


def test_predicted_taxonomy_conflict_and_safe_update(tmp_path):
    database = ComplexPortalDatabase(tmp_path / "db")
    good = fixture(tmp_path / "9606.tsv")
    first = database.install_from_file(good)
    second = database.install_from_file(good)
    assert second != first and database.active_snapshot() == second
    with pytest.raises(ValueError, match="Predicted"):
        database.install_from_file(fixture(tmp_path / "9606_predicted.tsv"))
    with pytest.raises(ValueError, match="Predicted"):
        database.install_from_file(fixture(tmp_path / "predicted_evidence.tsv", predicted=True))
    with pytest.raises(ValueError, match="non-human"):
        database.install_from_file(fixture(tmp_path / "mouse.tsv", tax="10090"))
    with pytest.raises(ValueError, match="Conflicting stoichiometry"):
        database.install_from_file(fixture(tmp_path / "conflict.tsv", conflict=True))
    with pytest.raises(ComplexPortalCancelled):
        database.install_from_file(good, cancel_requested=lambda: True)
    assert database.active_snapshot() == second and database.validate_snapshot()


def test_corruption_detected(tmp_path):
    source = fixture(tmp_path / "9606.tsv")
    database = ComplexPortalDatabase(tmp_path / "db")
    snapshot = database.install_from_file(source)
    raw = snapshot / "raw" / "9606.tsv"
    raw.write_text("corrupt", encoding="utf-8")
    with pytest.raises(ValueError, match="raw file hash"):
        database.validate_snapshot()
    second = database.install_from_file(source)
    sqlite_path = second / "complex_portal.sqlite"
    sqlite_path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="SQLite hash"):
        database.validate_snapshot()
    third = database.install_from_file(source)
    manifest_path = third / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scope"] = "predicted"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="scope"):
        database.validate_snapshot()
