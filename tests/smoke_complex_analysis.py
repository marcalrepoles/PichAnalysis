"""Controlled offline Python → R → Complex Portal scientific smoke."""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from pichanalysis.core.complex_analysis import ComplexParameters, run_complex_analysis
from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.r_runtime import RRuntime
from test_complex_portal_database import HEADERS


def build_fixture(root):
    source = root / "9606.tsv"
    def row(id, name, direct, expanded, evidence="ECO:0000353(manual assertion)"):
        result = dict.fromkeys(HEADERS, "-")
        result.update({"#Complex ac": id, "Recommended name": name, "Taxonomy identifier": "9606",
            "Identifiers (and stoichiometry) of molecules in complex": direct,
            "Expanded participant list": expanded, "Evidence Code": evidence,
            "Source": 'psi-mi:"MI:0469"(IntAct)'})
        return result
    rows = [
        row("CPX-1", "Alternative and nested", "P12345(4)|CPX-2(0)|CHEBI:123(1)|URS000001(0)", "P12345(4)|[Q11111,Q22222,Q33333](0)|Q44444(1)"),
        row("CPX-2", "Complete simple", "P12345(1)|Q44444(2)", "P12345(1)|Q44444(2)"),
        row("CPX-3", "Partial simple", "Q55555(1)|Q66666(1)", "Q55555(1)|Q66666(1)"),
        row("CPX-4", "Zero coverage", "Q77777(1)|Q88888(1)", "Q77777(1)|Q88888(1)"),
        row("CPX-5", "No protein", "CHEBI:124(1)", "-"),
        row("CPX-6", "Enriched", "P12345(1)|Q11111(1)|Q22222(1)|Q44444(1)|Q55555(1)", "P12345(1)|Q11111(1)|Q22222(1)|Q44444(1)|Q55555(1)"),
        row("CPX-7", "Non-significant", "P12345(1)|Q11111(1)|Q22222(1)|Q66666(1)|Q77777(1)|Q88888(1)", "P12345(1)|Q11111(1)|Q22222(1)|Q66666(1)|Q77777(1)|Q88888(1)"),
    ]
    with source.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    manager = DatabaseManager(root / "db")
    manager.complex_portal.install_from_file(source)
    manager.complex_portal.provider = None
    project = SimpleNamespace(root=root / "project", config={"organism_tax_id": "9606"})
    catalog = project.root / "mapping/tables/protein_catalog.csv"
    catalog.parent.mkdir(parents=True)
    accessions = ["P12345", "P12345-2", "Q11111", "Q22222", "Q44444", "Q55555", "Q66666", "Q77777", "Q88888", "O12345", "P12345;Q11111", "NONE"]
    pd.DataFrame([{"source_row": str(i), "original_id": f"row{i}", "uniprot_accession": value, "gene_symbol": f"Gene{i}"} for i,value in enumerate(accessions,1)]).to_csv(catalog,index=False)
    return project, manager


def main():
    root = Path(tempfile.mkdtemp(prefix="pichanalysis-complex-analysis-"))
    project, manager = build_fixture(root)
    outputs = run_complex_analysis(project, manager, RRuntime(), run_id="scientific_smoke",
        parameters=ComplexParameters(target_selection="Manual selection", manual_rows=(1,2,3,4,5,6)), timeout=180)
    tables = outputs["tables"]
    coverage = tables["complex_coverage"].set_index("complex_id")
    frequency = tables["complex_frequency"].set_index("complex_id")
    groups = tables["component_groups"]
    alternative = groups[(groups.complex_id == "CPX-1") & (groups.component_type == "alternative")].iloc[0]
    assert alternative.group_covered and alternative.option_count == 3
    assert alternative.detected_uniprot_accessions == "Q11111;Q22222"
    assert coverage.loc["CPX-1", "covered_component_groups"] == 3
    assert coverage.loc["CPX-1", "total_protein_component_groups"] == 3
    assert frequency.loc["CPX-1", "target_member_count"] == 4
    assert coverage.loc["CPX-2", "coverage_class"] == "complete_protein_component_coverage"
    assert coverage.loc["CPX-3", "coverage_class"] == "partial_protein_component_coverage"
    assert coverage.loc["CPX-4", "coverage_class"] == "no_detected_protein_components"
    assert coverage.loc["CPX-5", "coverage_class"] == "not_applicable_no_protein_components"
    assert coverage.loc["CPX-1", "has_nonprotein_participants"] and coverage.loc["CPX-1", "has_nested_complexes"]
    assert coverage.loc["CPX-1", "has_unknown_stoichiometry"]
    assert outputs["metadata"]["target_size"] == 5
    assert outputs["metadata"]["background_size"] == 9
    assert tables["mapping"].loc[tables["mapping"].source_row.isin((1,2)), "canonical_uniprot"].nunique() == 1
    assert len(tables["ambiguous"]) == 1 and len(tables["unmapped"]) == 1
    assert "CPX-6" in set(tables["complex_enrichment_significant"].complex_id)
    assert "CPX-7" not in set(tables["complex_enrichment_significant"].complex_id)
    assert "CPX-5" in set(tables["complex_enrichment_excluded"].complex_id)
    assert len(outputs["plots"]) == 12 and outputs["workbook"].is_file()
    required_sheets = {"Summary", "Protein mapping", "Ambiguous", "Unmapped", "Complex coverage", "Coverage classes", "Component groups", "Alternative groups", "Direct participants", "Frequency", "Enrichment", "Enrichment excluded", "Protein to complexes", "Complex to proteins", "Non-protein summary", "Nested complexes", "Stoichiometry"}
    assert required_sheets.issubset(set(pd.ExcelFile(outputs["workbook"]).sheet_names))
    assert (project.root / "scripts/runs/scientific_smoke_complexes/09_complex_analysis.R").is_file()
    assert (project.root / "scripts/runs/scientific_smoke_complexes/complex_analysis.R").is_file()
    assert (project.root / "scripts/runs/scientific_smoke_complexes/R_session_info.txt").is_file()
    assert manager.complex_portal.validate_snapshot()
    print({"temporary_path": str(root), "run_id": "scientific_smoke", "snapshot_id": outputs["metadata"]["snapshot_id"],
           "target_size": outputs["metadata"]["target_size"], "background_size": outputs["metadata"]["background_size"],
           "alternative_group": "[Q11111,Q22222,Q33333]", "detected_alternatives": "Q11111;Q22222",
           "coverage_contribution": 1, "frequency_enrichment_contribution": 2,
           "tested_complexes": len(tables["complex_enrichment_all"]),
           "significant_complexes": len(tables["complex_enrichment_significant"]), "offline": True})


if __name__ == "__main__":
    main()
