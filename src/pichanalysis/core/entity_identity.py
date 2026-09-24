"""Typed, conservative identities for navigation across frozen analysis runs.

An experimental feature is never equated with a protein or gene. No identifier
is inferred or globally canonicalized here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class IdentityType(StrEnum):
    EXPERIMENTAL_FEATURE = "experimental_feature"
    SOURCE_ROW = "source_row"
    UNIPROT = "uniprot"
    GENE_SYMBOL = "gene_symbol"
    NCBI_GENE = "ncbi_gene"
    ENSEMBL_GENE = "ensembl_gene"
    ENSEMBL_PROTEIN = "ensembl_protein"
    REFSEQ = "refseq"
    MODULE_SPECIFIC = "module_specific"


class MatchStatus(StrEnum):
    EXACT = "Exact"
    MAPPED = "Mapped"
    AMBIGUOUS = "Ambiguous"
    NOT_FOUND = "Not found"
    NO_COMPATIBLE_RUN = "No compatible run"
    INCOMPLETE_RUN = "Incomplete run"


@dataclass(frozen=True)
class InputLineage:
    project_id: str
    frozen_input_hash: str | None = None
    matrix_hash: str | None = None
    source: str = "Unknown"

    def compatible(self, other: "InputLineage") -> bool:
        if self.project_id != other.project_id:
            return False
        return bool(self.frozen_input_hash and other.frozen_input_hash
            and self.frozen_input_hash == other.frozen_input_hash)


@dataclass(frozen=True)
class EntityContext:
    source_module: str
    source_run_id: str
    source_result_type: str
    feature_id: str | None = None
    source_row: int | None = None
    original_identifier: str | None = None
    original_identifier_type: str | None = None
    uniprot_accessions: tuple[str, ...] = ()
    gene_symbols: tuple[str, ...] = ()
    ncbi_gene_ids: tuple[str, ...] = ()
    ensembl_gene_ids: tuple[str, ...] = ()
    ensembl_protein_ids: tuple[str, ...] = ()
    refseq_ids: tuple[str, ...] = ()
    module_specific_ids: tuple[tuple[str, str], ...] = ()
    protein_group_members: tuple[str, ...] = ()
    input_lineage: InputLineage | None = None
    mapping_provenance: str | None = None
    database_snapshot_provenance: str | None = None


@dataclass(frozen=True)
class CrossModuleMatch:
    target_module: str
    target_run_id: str
    target_entity: str
    status: MatchStatus
    basis: str
    source_identifier: str = ""
    target_identifier: str = ""
    lineage_compatible: bool = False
    database_snapshot: str = ""
    run_date: str = ""
    artifact_source: str = ""
    details: dict[str, str] = field(default_factory=dict)


def tokens(value: object) -> tuple[str, ...]:
    """Split only explicitly persisted multi-value fields; preserve isoforms."""
    import re
    if value is None:
        return ()
    text = str(value).strip()
    if not text or text.casefold() in {"na", "nan", "none", "null"}:
        return ()
    return tuple(dict.fromkeys(part for part in re.split(r"[;,|\s]+", text)
        if part and part.casefold() not in {"na", "nan", "null"}))
