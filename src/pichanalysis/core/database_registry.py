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

DATABASE_REGISTRY = {KEGG_HSA.database_id: KEGG_HSA}

