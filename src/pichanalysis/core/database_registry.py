from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class DatabaseState(StrEnum):
    NOT_INSTALLED = "Not installed"
    DOWNLOADING = "Downloading"
    INCOMPLETE = "Incomplete"
    READY = "Ready"
    ERROR = "Error"
    UPDATING = "Updating"
    CANCELLED = "Cancelled"


@dataclass(frozen=True)
class DatabaseSpec:
    database_id: str
    display_name: str
    organism_name: str
    organism_code: str
    taxonomy_id: str
    components: tuple[str, ...]


KEGG_HSA = DatabaseSpec(
    database_id="kegg",
    display_name="KEGG Homo sapiens",
    organism_name="Homo sapiens",
    organism_code="hsa",
    taxonomy_id="9606",
    components=("core", "entries", "kgml", "images"),
)
REACTOME_HUMAN=DatabaseSpec("reactome","Reactome — Homo sapiens","Homo sapiens","R-HSA-","9606",("core","diagrams"))
MITOCARTA_HUMAN=DatabaseSpec("mitocarta","MitoCarta3.0 — Homo sapiens","Homo sapiens","human","9606",("core",))
STRING_HUMAN=DatabaseSpec("string","STRING — Homo sapiens","Homo sapiens","human","9606",("functional","physical","proteins","aliases"))
COMPLEX_PORTAL_HUMAN=DatabaseSpec("complex_portal","Complex Portal — Homo sapiens","Homo sapiens","human","9606",("curated_complexes",))
GENE_ONTOLOGY_HUMAN=DatabaseSpec("gene_ontology","Gene Ontology — Homo sapiens","Homo sapiens","human","9606",("ontology","annotations"))
MTDNA_EVIDENCE_HUMAN=DatabaseSpec("mtdna_evidence","mtDNA Evidence — Homo sapiens","Homo sapiens","human","9606",("derived",))

DATABASE_REGISTRY = {KEGG_HSA.database_id: KEGG_HSA,REACTOME_HUMAN.database_id:REACTOME_HUMAN,MITOCARTA_HUMAN.database_id:MITOCARTA_HUMAN,STRING_HUMAN.database_id:STRING_HUMAN,COMPLEX_PORTAL_HUMAN.database_id:COMPLEX_PORTAL_HUMAN,MTDNA_EVIDENCE_HUMAN.database_id:MTDNA_EVIDENCE_HUMAN,GENE_ONTOLOGY_HUMAN.database_id:GENE_ONTOLOGY_HUMAN}
