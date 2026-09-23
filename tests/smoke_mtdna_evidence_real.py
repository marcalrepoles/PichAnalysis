"""Opt-in real-source mtDNA evidence smoke; never collected by normal pytest.

Downloads official MitoCarta, GO, and NCBI sources into an isolated temporary root.
"""
from __future__ import annotations

import argparse
import json
import tempfile
from collections import Counter
from pichanalysis.core.mtdna_evidence_database import GO_SEEDS, parse_ncbi_refseq, mitocarta_category
from pathlib import Path

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.databases.mitocarta import FILES, MitoCartaProvider
from pichanalysis.core.mtdna_evidence_database import MtdnaEvidenceDatabase


def main():
    temporary=Path(tempfile.mkdtemp(prefix="pichanalysis-mtdna-evidence-real-"))
    manager=DatabaseManager(temporary/"databases")
    mito=manager.mitocarta
    staging=mito.create_staging_snapshot()
    provider=MitoCartaProvider()
    files=[]
    for name in FILES:
        print(f"Downloading official MitoCarta: {name}", flush=True)
        files.append(provider.download(name,staging/"raw"/name))
    mito_snapshot=mito.finalize_snapshot(staging,files)
    print(f"MitoCarta snapshot Ready: {mito_snapshot}", flush=True)
    go_snapshot=manager.gene_ontology.download_and_build(
        progress=lambda payload: print(payload.get("message",""), flush=True))
    print(f"Gene Ontology snapshot Ready: {go_snapshot}", flush=True)
    derived=manager.mtdna_evidence
    snapshot=derived.build(mitocarta=mito,go_snapshot=go_snapshot)
    assert mito.is_ready() and manager.gene_ontology.is_ready()
    assert derived.validate_snapshot(snapshot)
    manifest=derived.manifest(snapshot)
    agreement=derived.source_agreement()
    mito_evidence=derived._query("SELECT evidence_category, COUNT(DISTINCT entity_key) AS genes FROM mitocarta_evidence GROUP BY evidence_category")
    go_evidence=derived._query("SELECT evidence_category, COUNT(DISTINCT entity_key) AS genes FROM go_evidence GROUP BY evidence_category")
    mito_categories=dict(zip(mito_evidence.evidence_category,mito_evidence.genes.astype(int)))
    go_categories=dict(zip(go_evidence.evidence_category,go_evidence.genes.astype(int)))
    if not any(name in mito_categories for name in ("mtDNA_replication","mitochondrial_nucleoid","mtDNA_repair")):
        raise AssertionError("Real MitoCarta source lacks expected mtDNA maintenance descendants.")
    if not any(name in go_categories for name in ("mitochondrial_nucleoid","mtDNA_metabolism","mtDNA_replication","mtDNA_repair","mitochondrial_transcription")):
        raise AssertionError("Real GO source lacks expected mitochondrial-specific annotations.")
    def example(frame):
        if frame.empty:return None
        key=sorted(frame.entity_key.astype(str))[0]
        row=frame.loc[frame.entity_key.eq(key)].iloc[0]
        evidence=derived.get_evidence_for_entity(key)
        entity=derived.lookup_gene(key).iloc[0]
        return {"entity_key":key,"gene":entity.gene_symbol,
            "mitocarta":bool(row.mitocarta_evidence),"go":bool(row.go_evidence),
            "mtdna_encoded":bool(row.mtdna_encoded),"source_count":int(row.source_count),
            "categories":sorted(evidence.evidence_category.unique().tolist()),
            "evidence_record_count":len(evidence),
            "evidence":evidence[["source","evidence_category","source_identifier","evidence_code"]].to_dict("records")}
    both=agreement[agreement.mitocarta_evidence.astype(bool)&agreement.go_evidence.astype(bool)]
    mito_only=agreement[agreement.mitocarta_evidence.astype(bool)&~agreement.go_evidence.astype(bool)]
    go_only=agreement[agreement.go_evidence.astype(bool)&~agreement.mitocarta_evidence.astype(bool)]
    encoded=agreement[agreement.mtdna_encoded.astype(bool)]
    report={"temporary_path":str(temporary),"mtdna_evidence_snapshot_id":snapshot.name,
        "mitocarta_snapshot_id":manifest["mitocarta_snapshot_id"],"mitocarta_hashes":manifest["mitocarta_source_hashes"],
        "go_snapshot_id":manifest["go_snapshot_id"],"go_hashes":manifest["go_source_hashes"],
        "ncbi_accession":manifest["ncbi_accession"],"ncbi_source_url":manifest["ncbi_source_url"],
        "ncbi_raw_size":manifest["ncbi_raw_size"],"ncbi_raw_sha256":manifest["ncbi_raw_sha256"],
        "entity_count":manifest["counts"]["entities"],"evidence_record_count":manifest["counts"]["evidence_records"],
        "mitocarta_supported_count":manifest["counts"]["mitocarta_supported_genes"],
        "go_supported_count":manifest["counts"]["go_supported_genes"],
        "mtdna_encoded_count":manifest["counts"]["mtdna_encoded_genes"],
        "category_counts":manifest["category_counts"],
        "mitocarta_category_counts":mito_categories,"go_category_counts":go_categories,
        "extracted_ncbi_protein_cds":manifest["counts"]["mtdna_encoded"],"sqlite_size":(snapshot/"mtdna_evidence.sqlite").stat().st_size,
        "sqlite_sha256":manifest["sqlite_sha256"],"examples":{"mitocarta_and_go":example(both),
            "mitocarta_only":example(mito_only),"go_only":example(go_only),"mtdna_encoded":example(encoded)}}
    derived.provider=None
    assert derived.validate_snapshot(snapshot)
    if not agreement.empty:
        key=agreement.entity_key.iloc[0]
        assert not derived.get_evidence_for_entity(key).empty
        assert not derived.lookup_gene(key).empty
        assert not derived.source_agreement(key).empty
    go_manifest=manager.gene_ontology.manifest(go_snapshot)
    report["go_ontology_version"]=go_manifest["ontology_version"]
    report["go_annotation_version"]=go_manifest["annotation_version"]
    report["go_raw_hashes"]={name:value["sha256"] for name,value in go_manifest["sources"].items()}
    report["go_counts"]=go_manifest["counts"]
    report["go_excluded_taxa"]=go_manifest["excluded_taxa"]
    report["mitocarta_category_details"]=derived._query(
        "SELECT source_path, evidence_category, COUNT(DISTINCT entity_key) AS entity_count, "
        "COUNT(*) AS evidence_record_count FROM mitocarta_evidence "
        "GROUP BY source_path, evidence_category ORDER BY source_path").to_dict("records")
    report["go_category_details"]=derived._query(
        "SELECT evidence_category, COUNT(DISTINCT entity_key) AS entity_count, "
        "COUNT(*) AS evidence_record_count FROM go_evidence GROUP BY evidence_category "
        "ORDER BY evidence_category").to_dict("records")
    report["category_membership_count"]=manifest["counts"]["category_membership"]
    report["multiple_source_entities"]=int((agreement.source_count.astype(int)>1).sum())
    report["symbol_fallback_entities"]=int(agreement.entity_key.str.startswith("SYMBOL:").sum())
    products=parse_ncbi_refseq((snapshot/"raw"/"NC_012920.1.xml").read_bytes())
    report["ncbi_cds_count"]=len(products)
    report["ncbi_gene_ids_found"]=sum(bool(row["ncbi_gene_id"]) for row in products)
    report["ncbi_protein_ids_found"]=sum(bool(row["protein_id"]) for row in products)
    report["ncbi_protein_products"]=sum(bool(row["product"]) for row in products)
    hierarchy=mito.table_path("hierarchy")
    import pandas as pd
    paths=pd.read_csv(hierarchy,dtype=str).fillna("")
    report["mitocarta_relevant_paths"]={row.pathway_full:mitocarta_category(row.pathway_full)
        for row in paths.itertuples(index=False) if mitocarta_category(row.pathway_full)}
    assert all("translation" not in path.casefold() and "ribosome" not in path.casefold()
        for path in report["mitocarta_relevant_paths"])
    report["go_seed_descendants"]={seed:int(len(manager.gene_ontology.get_descendants(seed,snapshot=go_snapshot))-1)
        for seed in GO_SEEDS}
    report["go_seed_positive_annotations"]={}
    for seed in GO_SEEDS:
        ids=manager.gene_ontology.get_descendants(seed,snapshot=go_snapshot).go_id.tolist()
        annotations=manager.gene_ontology.get_annotations_for_terms(ids,snapshot=go_snapshot)
        report["go_seed_positive_annotations"][seed]=int((~annotations.qualifier.str.contains(r"(?:^|\|)NOT(?:$|\|)",regex=True)).sum())
    for seed in GO_SEEDS:
        assert not manager.gene_ontology.get_term(seed,snapshot=go_snapshot).empty
    assert not both.empty and all(key.startswith("NCBI:") for key in both.entity_key)
    report["selected_go_not_exclusions"]=manifest["counts"]["go_exclusions"]
    with __import__("sqlite3").connect(snapshot/"mtdna_evidence.sqlite") as connection:
        accession=connection.execute("SELECT uniprot FROM evidence_records WHERE uniprot!='' ORDER BY uniprot LIMIT 1").fetchone()[0].split(";")[0]
    lookup=derived.lookup_uniprot(accession)
    assert lookup["status"] in ("unique","ambiguous") and not lookup["entities"].empty
    report["uniprot_lookup"]={"accession":accession,"status":lookup["status"],
        "entities":lookup["entities"].entity_key.tolist()}
    derived=MtdnaEvidenceDatabase(manager.root)
    derived.provider=None
    key=sorted(agreement.entity_key.astype(str))[0]
    assert not derived.lookup_gene(key).empty
    assert not derived.get_evidence_for_entity(key).empty
    assert not derived.get_categories_for_entity(key).empty
    assert not derived.get_entities_for_category(manifest["category_counts"].keys().__iter__().__next__()).empty
    derived.get_mitocarta_evidence(key)
    derived.get_go_evidence(key)
    derived.is_mtdna_encoded(key)
    assert not derived.search_entities(key.split(":")[-1]).empty
    assert not derived.source_agreement(key).empty
    assert derived.validate_snapshot(snapshot)
    report["offline_queries_verified"]=True
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
