"""Offline Python → R mtDNA Evidence smoke; no network or user database state."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pandas as pd

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.mtdna_analysis import MtdnaParameters, run_mtdna_analysis
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime
from test_mtdna_evidence_database import _xml, sources


def build_fixture(root):
    root = Path(root)
    project = create_project(root / "projects", "mtDNA smoke")
    set_organism(project, "Homo sapiens", "9606")
    manager = DatabaseManager(root / "databases")
    mito, go = sources(manager.root, "fixture-A")
    annotations = pd.read_csv(go / "tables/go_annotations.csv", dtype=str).fillna("")
    extra = pd.DataFrame([
        dict(ncbi_gene_id="2", gene_symbol="GENEB", uniprot="P22222", go_id="GO:0042645",
             evidence_code="IDA", reference="PMID:2", qualifier="", assigned_by="UniProt"),
        dict(ncbi_gene_id="3", gene_symbol="GENEC", uniprot="P33333", go_id="GO:0042645",
             evidence_code="IDA", reference="PMID:3", qualifier="", assigned_by="UniProt"),
    ])
    pd.concat([annotations, extra], ignore_index=True).to_csv(go / "tables/go_annotations.csv", index=False)
    ncbi = root / "ncbi.xml"
    ncbi.write_bytes(_xml(gene="MT-ATP6", gene_id="100"))
    manager.mtdna_evidence.build(mitocarta=mito, go_snapshot=go, ncbi_xml=ncbi)
    catalog = []
    for row in range(1, 21):
        catalog.append(dict(source_row=row, original_id=f"row{row}", mapping_status="mapped_unique",
            ncbi_gene_id=str(row), gene_symbol=f"GENE{row}" if row > 3 else f"GENE{chr(64+row)}",
            uniprot_accession=f"P{row:05d}"))
    catalog.append(dict(source_row=21, original_id="duplicate", mapping_status="mapped_unique",
        ncbi_gene_id="1", gene_symbol="GENEA", uniprot_accession="P11111-2"))
    catalog.append(dict(source_row=22, original_id="ambiguous", mapping_status="ambiguous",
        ncbi_gene_id="1;2", gene_symbol="GENEA;GENEB", uniprot_accession="P11111;P22222"))
    catalog.append(dict(source_row=23, original_id="mtDNA only", mapping_status="mapped_unique",
        ncbi_gene_id="100", gene_symbol="MT-ATP6", uniprot_accession=""))
    path = project.root / "mapping/tables/protein_catalog.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(catalog).to_csv(path, index=False)
    presence = project.root / "analyses/presence_absence/tables/classification.csv"
    presence.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"source_row": list(range(1, 24)),
        "classification": ["A-specific" if i <= 3 else "B-specific" if i <= 6 else "Shared"
                           for i in range(1, 24)]}).to_csv(presence, index=False)
    return project, manager


def main():
    temporary = Path(tempfile.mkdtemp(prefix="pichanalysis-mtdna-analysis-"))
    project, manager = build_fixture(temporary)
    runtime = RRuntime()
    outputs = run_mtdna_analysis(project, manager, runtime, run_id="offline-smoke",
        parameters=MtdnaParameters(target_selection="A-specific", minimum_overlap=2))
    assert outputs["metadata"]["network_access"] is False
    assert outputs["tables"]["enrichment_all"].shape[0] >= 1
    assert outputs["tables"]["enrichment_significant"].shape[0] >= 1
    assert outputs["tables"]["entity_evidence"].shape[0] >= 20
    assert outputs["workbook"].is_file()
    assert len(outputs["plots"]) >= 5
    print(json.dumps({"temporary_root": str(temporary), "run_id": "offline-smoke",
        "snapshot_id": outputs["metadata"]["snapshot_id"],
        "target_size": outputs["metadata"]["target_size"],
        "background_size": outputs["metadata"]["background_size"],
        "significant_categories": len(outputs["tables"]["enrichment_significant"]),
        "plots": len(outputs["plots"])}))


if __name__ == "__main__": main()
