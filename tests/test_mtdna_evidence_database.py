"""Offline source-specific mtDNA evidence and immutable snapshot regression."""
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.mtdna_evidence_database import (
    MtdnaBuildCancelled, MtdnaEvidenceDatabase, MtdnaEvidenceError,
    MissingMtdnaDependency, mitocarta_category, parse_ncbi_refseq,
)


def _xml(gene="GENEA",gene_id="1"):
    return f'''<GBSet><GBSeq><GBSeq_accession-version>NC_012920.1</GBSeq_accession-version>
<GBSeq_organism>Homo sapiens</GBSeq_organism><GBSeq_feature-table>
<GBFeature><GBFeature_key>source</GBFeature_key><GBFeature_quals><GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>taxon:9606</GBQualifier_value></GBQualifier></GBFeature_quals></GBFeature>
<GBFeature><GBFeature_key>CDS</GBFeature_key><GBFeature_location>1..100</GBFeature_location><GBFeature_quals>
<GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>{gene}</GBQualifier_value></GBQualifier>
<GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>GeneID:{gene_id}</GBQualifier_value></GBQualifier>
<GBQualifier><GBQualifier_name>product</GBQualifier_name><GBQualifier_value>protein A</GBQualifier_value></GBQualifier>
<GBQualifier><GBQualifier_name>protein_id</GBQualifier_name><GBQualifier_value>YP_000001.1</GBQualifier_value></GBQualifier>
</GBFeature_quals></GBFeature></GBSeq_feature-table></GBSeq></GBSet>'''.encode()


class LocalMito:
    def __init__(self,snapshot):self.snapshot=snapshot
    def active_snapshot(self):return self.snapshot
    def manifest(self,snapshot=None):return json.loads((self.snapshot/"manifest.json").read_text())


def sources(root,tag="A"):
    mito=root/"mitocarta"/"human"/"snapshots"/tag
    go=root/"go"/"human"/"snapshots"/tag
    (mito/"tables").mkdir(parents=True);(mito/"raw").mkdir()
    (go/"tables").mkdir(parents=True)
    raw=mito/"raw"/"Human.MitoCarta3.0.xls";raw.write_bytes(b"local-mitocarta-fixture")
    (mito/"manifest.json").write_text(json.dumps({"status":"Ready","files":[{"filename":raw.name,"sha256":hashlib.sha256(raw.read_bytes()).hexdigest()}]}))
    (go/"manifest.json").write_text(json.dumps({"status":"Ready","snapshot_id":tag,"tax_id":"9606"}))
    genes=pd.DataFrame([
        {"gene_symbol":"GENEA","ncbi_gene_id":"1","uniprot_accession":"P11111","is_mitocarta":True,"evidence":"curated","sub_compartment_raw":"matrix"},
        {"gene_symbol":"GENEB","ncbi_gene_id":"2","uniprot_accession":"P22222","is_mitocarta":True,"evidence":"curated","sub_compartment_raw":"matrix"},
        {"gene_symbol":"GENEC","ncbi_gene_id":"3","uniprot_accession":"P33333","is_mitocarta":True,"evidence":"curated","sub_compartment_raw":"matrix"},
        {"gene_symbol":"GRANULE","ncbi_gene_id":"4","uniprot_accession":"P44444","is_mitocarta":True,"evidence":"curated","sub_compartment_raw":"matrix"},
    ])
    genes.to_csv(mito/"tables/mitocarta_genes.csv",index=False)
    paths={"rep":"Mitochondrial central dogma > mtDNA maintenance > mtDNA replication",
        "nuc":"Mitochondrial central dogma > mtDNA maintenance > mtDNA nucleoid",
        "tx":"Mitochondrial central dogma > mtRNA metabolism > Transcription",
        "granule":"Mitochondrial central dogma > mtRNA metabolism > mtRNA granules"}
    pd.DataFrame([{"mitopathway":key,"pathway_full":value} for key,value in paths.items()]).to_csv(mito/"tables/mitopathway_hierarchy.csv",index=False)
    pd.DataFrame([{"mitopathway":"nuc","gene_symbol":"GENEA"},{"mitopathway":"rep","gene_symbol":"GENEB"},
        {"mitopathway":"tx","gene_symbol":"GENEC"},{"mitopathway":"granule","gene_symbol":"GRANULE"}]).to_csv(mito/"tables/mitopathway_membership.csv",index=False)
    terms={"GO:0042645":"mitochondrial nucleoid","GO:0006264":"mitochondrial DNA replication",
        "GO:0043504":"mitochondrial DNA repair","GO:0006390":"mitochondrial transcription",
        "GO:0006392":"mitochondrial transcription initiation","GO:1901858":"regulation of mitochondrial DNA metabolic process",
        "GO:1234567":"specific mtDNA regulation","GO:0006281":"DNA repair","GO:0006260":"DNA replication"}
    pd.DataFrame([{"go_id":key,"name":value} for key,value in terms.items()]).to_csv(go/"tables/go_terms.csv",index=False)
    pd.DataFrame([{"child_go_id":"GO:0006392","parent_go_id":"GO:0006390","relation":"part_of"},
        {"child_go_id":"GO:1234567","parent_go_id":"GO:1901858","relation":"is_a"}]).to_csv(go/"tables/go_relations.csv",index=False)
    base={"evidence_code":"IDA","reference":"PMID:1","qualifier":"","assigned_by":"UniProt"}
    rows=[]
    for gene,gid,uni,goid in (("GENEA","1","P11111","GO:0042645"),("GENEB","2","P22222","GO:0006264"),
        ("GENEC","3","P33333","GO:0043504"),("GENED","5","P55555","GO:0006392"),
        ("REG","6","P66666","GO:1234567"),("GENERIC","7","P77777","GO:0006281"),
        ("NOTGENE","8","P88888","GO:0042645"),("OTHER","9","P11111","GO:0006390")):
        rows.append(dict(base,ncbi_gene_id=gid,gene_symbol=gene,uniprot=uni,go_id=goid,
            qualifier="NOT|located_in" if gene=="NOTGENE" else ""))
    pd.DataFrame(rows).to_csv(go/"tables/go_annotations.csv",index=False)
    return LocalMito(mito),go


def test_parsing_and_scientific_exclusions():
    products=parse_ncbi_refseq(_xml())
    assert products==[{"gene_symbol":"GENEA","ncbi_gene_id":"1","protein_id":"YP_000001.1","product":"protein A","coordinates":"1..100"}]
    with pytest.raises(MtdnaEvidenceError):parse_ncbi_refseq(_xml().replace(b"taxon:9606",b"taxon:10090"))
    assert mitocarta_category("Mitochondrial central dogma > mtRNA metabolism > Transcription")=="mitochondrial_transcription"
    assert mitocarta_category("Mitochondrial central dogma > mtRNA metabolism > mtRNA granules") is None
    assert mitocarta_category("Mitochondrial central dogma > Translation") is None


def test_offline_build_queries_history_and_corruption(tmp_path):
    db=MtdnaEvidenceDatabase(tmp_path)
    with pytest.raises(MissingMtdnaDependency):db.build(ncbi_xml=tmp_path/"missing.xml")
    mito,go=sources(tmp_path,"A")
    ncbi=tmp_path/"ncbi.xml";ncbi.write_bytes(_xml())
    x=db.build(mitocarta=mito,go_snapshot=go,ncbi_xml=ncbi)
    assert db.is_ready() and db.validate_snapshot(x)
    m=db.manifest(x);assert m["mitocarta_snapshot_id"]==m["go_snapshot_id"]=="A"
    assert m["ncbi_accession"]=="NC_012920.1"
    evidence=db.get_evidence_for_entity("NCBI:1")
    assert len(evidence)==3 and set(evidence.source)=={"MitoCarta","Gene Ontology","NCBI mtDNA genome"}
    assert db.source_agreement("NCBI:1").source_count.iloc[0]==3
    assert db.is_mtdna_encoded("NCBI:1")
    assert set(db.get_categories_for_entity("NCBI:1").category)=={"mitochondrial_nucleoid","mtDNA_encoded"}
    assert not db.is_mtdna_encoded("NCBI:2")
    assert "mitochondrial_nucleoid" not in set(db.get_categories_for_entity("NCBI:2").category)
    assert db.get_go_evidence("NCBI:5").evidence_category.iloc[0]=="mitochondrial_transcription"
    assert db.get_go_evidence("NCBI:6").evidence_class.iloc[0]=="regulation"
    assert db.get_evidence_for_entity("NCBI:7").empty
    assert db.get_evidence_for_entity("NCBI:8").empty
    assert not db.get_entities_for_category("mtDNA_repair").empty
    assert len(db.lookup_gene("GENEA"))==1
    assert db.lookup_uniprot("P11111")["status"]=="ambiguous"
    assert not db.search_entities("nucleoid").empty
    assert pd.read_csv(x/"tables/go_exclusions.csv").qualifier.iloc[0]=="NOT|located_in"
    mito_b,go_b=sources(tmp_path,"B")
    y=db.build(mitocarta=mito_b,go_snapshot=go_b,ncbi_xml=ncbi)
    assert db.active_snapshot()==y and db.manifest(x)["go_snapshot_id"]=="A"
    assert db.lookup_gene("GENEA",snapshot=x).entity_key.iloc[0]=="NCBI:1"
    with pytest.raises(MtdnaBuildCancelled):db.build(mitocarta=mito_b,go_snapshot=go_b,ncbi_xml=ncbi,cancel_requested=lambda:True)
    assert db.active_snapshot()==y
    bad=tmp_path/"invalid.xml";bad.write_bytes(b"bad")
    with pytest.raises(MtdnaEvidenceError):db.build(mitocarta=mito_b,go_snapshot=go_b,ncbi_xml=bad)
    assert db.active_snapshot()==y
    go_manifest=go/"manifest.json"
    original_go=go_manifest.read_bytes()
    go_manifest.write_bytes(b"{invalid")
    with pytest.raises(MtdnaEvidenceError):db.validate_snapshot(x)
    go_manifest.write_bytes(original_go)
    sqlite=x/"mtdna_evidence.sqlite"
    original_sqlite=sqlite.read_bytes()
    sqlite.write_bytes(b"corrupt")
    with pytest.raises(MtdnaEvidenceError):db.validate_snapshot(x)
    sqlite.write_bytes(original_sqlite)
    (x/"raw/NC_012920.1.xml").write_bytes(b"corrupt")
    with pytest.raises(MtdnaEvidenceError):db.validate_snapshot(x)
