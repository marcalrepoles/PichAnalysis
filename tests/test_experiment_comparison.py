from types import SimpleNamespace

import pandas as pd
import pytest

from pichanalysis.core.experiment_comparison import (ComparisonError, ComparisonSpec,
    TableSpec, compare, compare_differential, load_run)


def project(tmp_path):
    return SimpleNamespace(root=tmp_path)


def fixture_file(tmp_path, name, values):
    path = tmp_path / name
    pd.DataFrame({"accession": values}).to_csv(path, index=False)
    return path


def test_frozen_sets_duplicates_isoforms_and_source_mutation(tmp_path):
    a = fixture_file(tmp_path, "a.csv", ["P12345", "P12345", "P12345-2", "Q11111"])
    b = fixture_file(tmp_path, "b.csv", ["P12345", "B22222"])
    spec = ComparisonSpec(TableSpec(str(a), "accession", "uniprot"),
        TableSpec(str(b), "accession", "uniprot"), "uniprot")
    root = compare(project(tmp_path), spec, run_id="A")
    loaded = load_run(project(tmp_path), "A")
    assert loaded["manifest"]["summary"]["A duplicate identifier rows"] == 1
    assert set(loaded["tables"]["a_only"].comparison_entity) == {"P12345-2", "Q11111"}
    assert set(loaded["tables"]["b_only"].comparison_entity) == {"B22222"}
    assert set(loaded["tables"]["shared"].comparison_entity) == {"P12345"}
    assert loaded["tables"]["shared"].iloc[0].A_source_rows == "1;2"
    a.unlink()
    b.write_text("accession\nCHANGED\n", encoding="utf-8")
    assert set(load_run(project(tmp_path), "A")["tables"]["shared"].comparison_entity) == {"P12345"}
    assert (root / "experiment_comparison.xlsx").is_file()


def test_mapping_ambiguity_and_unmapped(tmp_path):
    a = fixture_file(tmp_path, "a.csv", ["P12345", "P12346", "P12347"])
    b = tmp_path / "b.csv"
    pd.DataFrame({"gene": ["GENE1"]}).to_csv(b, index=False)
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12346"],
        "gene_symbol": ["GENE1", "GENE2", "GENE3"]}).to_csv(mapping, index=False)
    spec = ComparisonSpec(TableSpec(str(a), "accession", "uniprot"),
        TableSpec(str(b), "gene", "gene_symbol"), "gene_symbol", str(mapping))
    loaded = load_run(project(tmp_path), compare(project(tmp_path), spec).name)
    assert set(loaded["tables"]["shared"].comparison_entity) == {"GENE1"}
    assert len(loaded["tables"]["ambiguous"]) == 1
    assert len(loaded["tables"]["unmapped_a"]) == 1


def test_missing_mapping_is_visible(tmp_path):
    a = fixture_file(tmp_path, "a.csv", ["P12345"])
    b = fixture_file(tmp_path, "b.csv", ["GENE1"])
    spec = ComparisonSpec(TableSpec(str(a), "accession", "uniprot"),
        TableSpec(str(b), "accession", "gene_symbol"), "gene_symbol")
    loaded = load_run(project(tmp_path), compare(project(tmp_path), spec).name)
    assert loaded["tables"]["unmapped_a"].iloc[0].status == "Mapping required"


def test_xlsx_requires_explicit_sheet(tmp_path):
    a = tmp_path / "multi.xlsx"
    with pd.ExcelWriter(a) as writer:
        pd.DataFrame({"accession": ["P12345"]}).to_excel(writer, sheet_name="One", index=False)
        pd.DataFrame({"accession": ["Q11111"]}).to_excel(writer, sheet_name="Two", index=False)
    b = fixture_file(tmp_path, "b.csv", ["P12345"])
    spec = ComparisonSpec(TableSpec(str(a), "accession", "uniprot"),
        TableSpec(str(b), "accession", "uniprot"), "uniprot")
    with pytest.raises(Exception, match="worksheet"):
        compare(project(tmp_path), spec)
    selected = ComparisonSpec(TableSpec(str(a), "accession", "uniprot", sheet="Two"),
        spec.b, "uniprot")
    loaded = load_run(project(tmp_path), compare(project(tmp_path), selected).name)
    assert loaded["config"]["a"]["sheet"] == "Two"


def test_differential_confirmation_controls_classes():
    a = pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12347"],
        "combined_significant": ["TRUE", "TRUE", "FALSE"],
        "effect_direction": ["higher_in_condition_A", "higher_in_condition_B", "higher_in_condition_A"],
        "effect": ["1", "-1", "0.2"]})
    b = pd.DataFrame({"uniprot_accession": ["P12345", "P12346", "P12347"],
        "combined_significant": ["TRUE", "TRUE", "FALSE"],
        "effect_direction": ["higher_in_condition_A", "higher_in_condition_A", "higher_in_condition_B"],
        "effect": ["2", "1", "-0.1"]})
    entities = {"P12345", "P12346", "P12347"}
    unconfirmed, _ = compare_differential(a, b, entities, False)
    assert set(unconfirmed.comparison_class) == {"Not comparable"}
    confirmed, _ = compare_differential(a, b, entities, True)
    assert confirmed.set_index("comparison_entity").loc["P12345", "comparison_class"] == "Significant up in both"
    assert confirmed.set_index("comparison_entity").loc["P12346", "comparison_class"] == "Significant opposite direction"
    assert confirmed.set_index("comparison_entity").loc["P12347", "comparison_class"] == "Not significant in either"
