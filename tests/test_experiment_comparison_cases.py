from types import SimpleNamespace

import pandas as pd
import pytest

from pichanalysis.core.experiment_comparison import (ComparisonError, ComparisonSpec,
    TableSpec, compare, compare_differential, load_run, suggestions)


def write(tmp_path, name, column, values):
    path = tmp_path / name
    pd.DataFrame({column: values, "description": [f"row {i}" for i in range(len(values))]}).to_csv(
        path, sep="\t" if name.endswith("tsv") else ",", index=False)
    return path


def run(tmp_path, a, b, kind_a="uniprot", kind_b="uniprot", target="uniprot", mapping=None, run_id=None):
    project = SimpleNamespace(root=tmp_path)
    spec = ComparisonSpec(TableSpec(str(a), "id", kind_a), TableSpec(str(b), "id", kind_b),
        target, str(mapping) if mapping else None)
    return load_run(project, compare(project, spec, run_id=run_id).name)


def test_tsv_csv_and_exact_uniprot(tmp_path):
    a, b = write(tmp_path, "a.tsv", "id", ["P12345", "Q11111"]), write(
        tmp_path, "b.csv", "id", ["P12345", "B22222"])
    result = run(tmp_path, a, b)
    assert result["manifest"]["summary"]["Shared"] == 1
    assert result["tables"]["shared"].iloc[0].A_mapping_status == "Exact"


@pytest.mark.parametrize("identifier_type,comparison_type,values", [
    ("gene_symbol", "gene_symbol", ["GeneA", "GENEB"]),
    ("entrez", "ncbi_gene", ["123", "456"]),
    ("uniprot", "uniprot", ["P12345", "P12345-2"]),
])
def test_exact_types_and_isoforms(tmp_path, identifier_type, comparison_type, values):
    a, b = write(tmp_path, "a.csv", "id", values), write(tmp_path, "b.csv", "id", [values[0]])
    result = run(tmp_path, a, b, identifier_type, identifier_type, comparison_type)
    assert result["manifest"]["summary"]["Shared"] == 1
    assert result["manifest"]["summary"]["A only"] == 1


def test_inference_is_suggestion_and_override_is_explicit(tmp_path):
    frame = pd.DataFrame({"UniProt": ["P12345", "Q11111"], "Intensity A1": [10, 20]})
    inferred = suggestions(frame)
    assert inferred["UniProt"]["identifier_type"] == "uniprot"
    assert inferred["Intensity A1"]["role"] == "quantification"
    a, b = write(tmp_path, "a.csv", "id", ["P12345"]), write(tmp_path, "b.csv", "id", ["P12345"])
    result = run(tmp_path, a, b, "original", "original", "original")
    assert result["manifest"]["summary"]["Shared"] == 1


def test_missing_column_and_empty_table_rejected(tmp_path):
    a, b = write(tmp_path, "a.csv", "id", ["P12345"]), write(tmp_path, "b.csv", "id", ["P12345"])
    spec = ComparisonSpec(TableSpec(str(a), "wrong", "uniprot"),
        TableSpec(str(b), "id", "uniprot"), "uniprot")
    with pytest.raises(ComparisonError, match="column"):
        compare(SimpleNamespace(root=tmp_path), spec)
    a.write_text("id,description\n", encoding="utf-8")
    with pytest.raises(Exception, match="no usable data"):
        run(tmp_path, a, b)


def test_mapping_provenance_ambiguous_gene_and_hash(tmp_path):
    a = write(tmp_path, "a.csv", "id", ["GENE1", "GENE2"])
    b = write(tmp_path, "b.csv", "id", ["P12345"])
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame({"gene_symbol": ["GENE1", "GENE2", "GENE2"],
        "uniprot_accession": ["P12345", "P12346", "P12347"]}).to_csv(mapping, index=False)
    result = run(tmp_path, a, b, "gene_symbol", "uniprot", "uniprot", mapping)
    assert result["manifest"]["summary"]["Ambiguous"] == 1
    assert result["tables"]["shared"].iloc[0].A_mapping_status == "Mapped"
    assert result["tables"]["shared"].iloc[0].mapping_provenance == str(mapping)
    frozen = result["run_root"] / "inputs/mapping/mapping.csv"
    frozen.write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ComparisonError, match="Mapping catalog integrity"):
        load_run(SimpleNamespace(root=tmp_path), result["run_root"].name)


def test_group_membership_not_duplicate_rows(tmp_path):
    a, b = write(tmp_path, "a.csv", "id", ["P12345;P12346"]), write(
        tmp_path, "b.csv", "id", ["P12346"])
    result = run(tmp_path, a, b)
    assert result["manifest"]["summary"]["Table A original rows"] == 1
    assert result["manifest"]["summary"]["A unique comparison entities"] == 2
    assert "group membership" in result["tables"]["shared"].iloc[0].match_basis


def test_two_history_runs_do_not_borrow_after_project_or_source_mutation(tmp_path):
    a, b = write(tmp_path, "a.csv", "id", ["P12345"]), write(tmp_path, "b.csv", "id", ["P12345"])
    first = run(tmp_path, a, b, run_id="first")
    a = write(tmp_path, "a.csv", "id", ["Q11111"])
    b = write(tmp_path, "b.csv", "id", ["B22222"])
    second = run(tmp_path, a, b, run_id="second")
    project = SimpleNamespace(root=tmp_path, config={"input": {"processed_file": "changed.csv"}})
    assert len(load_run(project, "first")["tables"]["shared"]) == 1
    assert len(load_run(project, "second")["tables"]["shared"]) == 0
    a.unlink()
    b.unlink()
    assert len(load_run(project, "first")["tables"]["shared"]) == 1
    assert first["run_root"] != second["run_root"]


@pytest.mark.parametrize("a_dir,a_sig,b_dir,b_sig,expected", [
    ("higher_in_condition_A", "TRUE", "higher_in_condition_A", "TRUE", "Significant up in both"),
    ("higher_in_condition_B", "TRUE", "higher_in_condition_B", "TRUE", "Significant down in both"),
    ("higher_in_condition_A", "TRUE", "higher_in_condition_B", "TRUE", "Significant opposite direction"),
    ("higher_in_condition_B", "TRUE", "higher_in_condition_A", "TRUE", "Significant opposite direction"),
    ("higher_in_condition_A", "TRUE", "higher_in_condition_B", "FALSE", "Significant only in A"),
    ("higher_in_condition_A", "FALSE", "higher_in_condition_B", "TRUE", "Significant only in B"),
    ("higher_in_condition_A", "FALSE", "higher_in_condition_B", "FALSE", "Not significant in either"),
])
def test_frozen_differential_patterns(a_dir, a_sig, b_dir, b_sig, expected):
    def frame(direction, significant):
        return pd.DataFrame({"uniprot_accession": ["P12345"], "effect_direction": [direction],
            "combined_significant": [significant], "effect": ["1.2"],
            "P.Value": ["0.01"], "adj.P.Val": ["0.02"]})
    result, _ = compare_differential(frame(a_dir, a_sig), frame(b_dir, b_sig), {"P12345"}, True)
    assert result.iloc[0].comparison_class == expected
    assert result.iloc[0]["A_effect"] == "1.2"


def test_na_not_tested_and_unconfirmed_contrasts():
    a = pd.DataFrame({"uniprot_accession": ["P12345"], "effect_direction": ["higher_in_condition_A"],
        "combined_significant": ["TRUE"], "effect": ["1"]})
    empty = pd.DataFrame(columns=a.columns)
    result, _ = compare_differential(a, empty, {"P12345"}, True)
    assert result.iloc[0].B_direction == "Not tested"
    assert pd.isna(result.iloc[0]["B_effect"])
    unconfirmed, _ = compare_differential(a, a, {"P12345"}, False)
    assert unconfirmed.iloc[0].comparison_class == "Not comparable"
