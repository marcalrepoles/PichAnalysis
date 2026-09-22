"""Explicit online smoke of the official curated human ComplexTab; excluded from pytest."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from pichanalysis.core.complex_portal_database import ComplexPortalDatabase


def main():
    root = Path(tempfile.mkdtemp(prefix="pichanalysis-complex-portal-"))
    database = ComplexPortalDatabase(root / "databases")
    snapshot = database.download_and_build()
    manifest = database.manifest(snapshot)
    complex_id = database.search_complexes("", limit=1)[0]["complex_id"]
    complex_metadata = database.get_complex(complex_id)
    direct = database.get_direct_participants(complex_id)
    proteins = database.get_proteins_for_complex(complex_id)
    with_protein = database._query("SELECT uniprot_accession FROM protein_membership ORDER BY uniprot_accession LIMIT 1")[0]["uniprot_accession"]
    memberships = database.get_complexes_for_uniprot(with_protein)
    isoform = database.get_complexes_for_uniprot(with_protein + "-2")
    known = database._query("SELECT * FROM direct_participants WHERE stoichiometry_known='1' LIMIT 1")
    unknown = database._query("SELECT * FROM direct_participants WHERE stoichiometry_known='0' LIMIT 1")
    nested = database._query("SELECT parent_complex_id FROM nested_complexes LIMIT 1")
    nonprotein = database._query("SELECT complex_id FROM nonprotein_participants LIMIT 1")
    database.provider = None
    assert database.validate_snapshot()
    assert database.get_complex(complex_id) == complex_metadata
    assert database.search_complexes(complex_id)
    assert database.get_complexes_for_uniprot(with_protein) == memberships
    assert database.get_proteins_for_complex(complex_id) == proteins
    assert database.get_direct_participants(complex_id) == direct
    assert database.get_nonprotein_participants(complex_id) is not None
    assert database.get_nested_complexes(complex_id) is not None
    assert memberships and isoform and isoform[0]["lookup_accession"] == with_protein
    assert isoform[0]["input_accession"] == with_protein + "-2" and isoform[0]["isoform_normalized"]
    assert known and known[0]["stoichiometry"] not in (None, "")
    if unknown:
        assert unknown[0]["stoichiometry"] == "" and unknown[0]["stoichiometry_raw"] == "0"
    if nested:
        assert database.get_nested_complexes(nested[0]["parent_complex_id"])
    if nonprotein:
        assert database.get_nonprotein_participants(nonprotein[0]["complex_id"])
    report = {
        "temporary_path": str(root), "snapshot_id": snapshot.name,
        "complex_portal_release": manifest["complex_portal_release"],
        "source_url": manifest["source_url"],
        "last_modified": manifest["http_last_modified"], "etag": manifest["http_etag"],
        "raw_filename": manifest["raw_file"]["filename"],
        "raw_size": manifest["raw_file"]["size"],
        "raw_sha256": manifest["raw_file"]["sha256"],
        "detected_headers": manifest["detected_headers"],
        "complex_count": manifest["complex_count"],
        "unique_protein_count": manifest["unique_protein_count"],
        "protein_membership_count": manifest["protein_membership_count"],
        "direct_participant_count": manifest["direct_participant_count"],
        "nonprotein_participant_count": manifest["nonprotein_participant_count"],
        "nested_complex_count": manifest["nested_complex_count"],
        "known_stoichiometry_count": manifest["known_stoichiometry_count"],
        "unknown_stoichiometry_count": manifest["unknown_stoichiometry_count"],
        "sqlite_size": manifest["sqlite_size"], "sqlite_sha256": manifest["sqlite_sha256"],
        "complex_lookup": {"complex_id": complex_id, "recommended_name": complex_metadata["recommended_name"],
                           "direct_participants": len(direct), "expanded_proteins": len(proteins)},
        "protein_lookup": {"uniprot_accession": with_protein, "complex_count": len(memberships)},
        "known_stoichiometry_example": known[0] if known else None,
        "unknown_stoichiometry_example": unknown[0] if unknown else None,
        "nested_example": nested[0] if nested else "No nested complex participant found in current curated human snapshot.",
        "nonprotein_example": nonprotein[0] if nonprotein else "No non-protein participant found in current curated human snapshot.",
        "offline_queries_passed": True,
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
