import pandas as pd

from pichanalysis.core.identifier_detection import IdentifierType, detect_identifier


def detect(name, values):
    return detect_identifier(name, pd.Series(values))


def test_detect_simple_uniprot():
    result = detect("Protein IDs", ["P12345", "Q9H0H5", "Q99999"])
    assert result.detected_type == IdentifierType.UNIPROT
    assert result.confidence >= 0.9


def test_detect_a0a_uniprot():
    assert detect("Accession", ["A0A123ABC4"]).detected_type == IdentifierType.UNIPROT


def test_detect_multiple_uniprot_in_cell():
    result = detect("Protein groups", ["P12345;Q99999", "Q9H0H5 P12345"])
    assert result.detected_type == IdentifierType.UNIPROT
    assert result.contains_multiple


def test_detect_ensembl_gene_with_version():
    assert detect("Ensembl gene", ["ENSG00000141510.18"]).detected_type == IdentifierType.ENSEMBL_GENE


def test_detect_ensembl_protein():
    assert detect("Ensembl protein", ["ENSP00000269305"]).detected_type == IdentifierType.ENSEMBL_PROTEIN


def test_detect_refseq_protein():
    assert detect("RefSeq", ["NP_000537.3", "XP_011527831.1"]).detected_type == IdentifierType.REFSEQ_PROTEIN


def test_numeric_column_is_not_aggressively_entrez():
    result = detect("Measurement", [100, 200, 300])
    assert result.detected_type == IdentifierType.UNKNOWN
    assert result.confidence < 0.55


def test_possible_gene_symbols_are_conservative():
    result = detect("Gene names", ["TFAM", "POLG", "NONO", "SFPQ", "LIG3"])
    assert result.detected_type == IdentifierType.GENE_SYMBOL
    assert result.confidence < 0.9
    assert "confirmation" in result.reason


def test_empty_column_is_unknown():
    result = detect("Protein", [None, float("nan")])
    assert result.detected_type == IdentifierType.UNKNOWN
    assert result.tested_values == 0

