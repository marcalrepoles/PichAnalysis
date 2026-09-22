import csv

import pytest

from pichanalysis.core.complex_portal_database import ComplexPortalDatabase
from pichanalysis.core.databases.complex_portal import ComplexPortalCancelled
from test_complex_portal_database import fixture


def rewrite(path, transform):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers, rows = reader.fieldnames, list(reader)
    transform(rows)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_alternative_protein_set_and_duplicate_membership(tmp_path):
    source = fixture(tmp_path / "9606.tsv")
    def change(rows):
        rows[0]["Expanded participant list"] += "|[P0DP23,P0DP24,P0DP25](1)|P12345(1)"
    rewrite(source, change)
    database = ComplexPortalDatabase(tmp_path / "db")
    database.install_from_file(source)
    expanded = database._query("SELECT * FROM expanded_protein_components WHERE complex_id='CPX-1'")
    assert len(expanded) == 5
    alternative = [row for row in expanded if row["participant_type"] == "protein_set"]
    assert len(alternative) == 1
    assert alternative[0]["participant_accession_raw"] == "[P0DP23,P0DP24,P0DP25]"
    assert alternative[0]["uniprot_accession"] == ""
    assert len(database.get_proteins_for_complex("CPX-1")) == 3
    assert not database.get_complexes_for_uniprot("P0DP23")
    assert database.validate_snapshot()


def test_no_protein_complex_is_retained_and_late_cancel_preserves_active(tmp_path):
    source = fixture(tmp_path / "9606.tsv")
    database = ComplexPortalDatabase(tmp_path / "db")
    first = database.install_from_file(source)
    def change(rows):
        rows[2]["Identifiers (and stoichiometry) of molecules in complex"] = "CHEBI:123(1)"
        rows[2]["Expanded participant list"] = "-"
    rewrite(source, change)
    canceled = False
    def progress(payload):
        nonlocal canceled
        if payload["stage"] == "index":
            canceled = True
    with pytest.raises(ComplexPortalCancelled):
        database.install_from_file(source, progress=progress, cancel_requested=lambda: canceled)
    assert database.active_snapshot() == first
    second = database.install_from_file(source)
    assert second != first
    assert database.get_complex("CPX-3") is not None
    assert database.get_proteins_for_complex("CPX-3") == []
    assert database.get_nonprotein_participants("CPX-3")[0]["participant_type"] == "chemical"
