import pytest

from pichanalysis.core.complex_portal_database import ComplexPortalDatabase
from test_complex_portal_database import fixture
from test_complex_portal_edge_cases import rewrite


def test_three_valid_complex_memberships_and_ready_immutability(tmp_path):
    source = fixture(tmp_path / "9606.tsv")
    rewrite(source, lambda rows: rows[1].update({
        "Identifiers (and stoichiometry) of molecules in complex": "P12345(1)|Q11111(2)|Q22222(4)",
        "Expanded participant list": "P12345(1)|Q11111(2)|Q22222(4)",
    }))
    database = ComplexPortalDatabase(tmp_path / "db")
    snapshot = database.install_from_file(source)
    assert [row["complex_id"] for row in database.get_complexes_for_uniprot("P12345")] == ["CPX-1", "CPX-2", "CPX-3"]
    with pytest.raises(ValueError, match="immutable"):
        database.mark_downloading(snapshot)
    with pytest.raises(ValueError, match="immutable"):
        database.mark_error(snapshot, "accidental")
    with pytest.raises(ValueError, match="immutable"):
        database.mark_incomplete(snapshot)
    assert database.validate_snapshot()
