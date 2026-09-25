from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import Iterable

import pandas as pd


MAX_TEST_VALUES = 200


class IdentifierType(str, Enum):
    UNIPROT = "uniprot"
    GENE_SYMBOL = "gene_symbol"
    ENTREZ = "entrez"
    ENSEMBL_GENE = "ensembl_gene"
    ENSEMBL_PROTEIN = "ensembl_protein"
    REFSEQ_PROTEIN = "refseq_protein"
    UNKNOWN = "unknown"


IDENTIFIER_LABELS = {
    IdentifierType.UNIPROT: "UniProt accession",
    IdentifierType.GENE_SYMBOL: "Gene symbol",
    IdentifierType.ENTREZ: "NCBI Gene ID / Entrez",
    IdentifierType.ENSEMBL_GENE: "Ensembl Gene ID",
    IdentifierType.ENSEMBL_PROTEIN: "Ensembl Protein ID",
    IdentifierType.REFSEQ_PROTEIN: "RefSeq Protein",
    IdentifierType.UNKNOWN: "Unknown",
}


@dataclass(frozen=True)
class DetectionResult:
    detected_type: IdentifierType
    confidence: float
    matched_values: int
    tested_values: int
    reason: str
    contains_multiple: bool = False

    @property
    def label(self) -> str:
        return IDENTIFIER_LABELS[self.detected_type]


_UNIPROT = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}){1,2}[0-9]|A0A[A-Z0-9]{7})$",
    re.IGNORECASE,
)
_ENSEMBL_GENE = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)
_ENSEMBL_PROTEIN = re.compile(r"^ENSP\d+(?:\.\d+)?$", re.IGNORECASE)
_REFSEQ = re.compile(r"^(?:NP|XP|YP|WP)_\d+(?:\.\d+)?$", re.IGNORECASE)
_GENE_SYMBOL = re.compile(r"^[A-Z][A-Z0-9-]{1,14}$")
_INTEGER = re.compile(r"^[1-9]\d*$")


def _values(series: pd.Series) -> list[str]:
    values: list[str] = []
    for value in series.dropna():
        text = str(value).strip()
        if text:
            values.append(text)
        if len(values) >= MAX_TEST_VALUES:
            break
    return values


def _tokens(value: str) -> list[str]:
    return [token for token in re.split(r"[;\s]+", value.strip()) if token]


def _ratio(values: Iterable[str], pattern: re.Pattern[str], split: bool = False) -> tuple[int, int, bool]:
    matched = 0
    tested = 0
    multiple = False
    for value in values:
        tokens = _tokens(value) if split else [value]
        multiple = multiple or len(tokens) > 1
        tested += 1
        if tokens and all(pattern.fullmatch(token) for token in tokens):
            matched += 1
    return matched, tested, multiple


def detect_identifier(column_name: str, series: pd.Series) -> DetectionResult:
    values = _values(series)
    if not values:
        return DetectionResult(IdentifierType.UNKNOWN, 0.0, 0, 0, "Column has no values to evaluate.")

    header = column_name.casefold()
    candidates: list[tuple[IdentifierType, re.Pattern[str], bool, tuple[str, ...]]] = [
        (IdentifierType.UNIPROT, _UNIPROT, True, ("uniprot", "protein id", "protein group", "accession")),
        (IdentifierType.ENSEMBL_GENE, _ENSEMBL_GENE, False, ("ensembl gene", "ensg")),
        (IdentifierType.ENSEMBL_PROTEIN, _ENSEMBL_PROTEIN, False, ("ensembl protein", "ensp")),
        (IdentifierType.REFSEQ_PROTEIN, _REFSEQ, False, ("refseq", "protein accession")),
    ]
    scored: list[DetectionResult] = []
    for kind, pattern, split, hints in candidates:
        matched, tested, multiple = _ratio(values, pattern, split)
        ratio = matched / tested
        header_bonus = 0.08 if any(hint in header for hint in hints) else 0.0
        confidence = min(1.0, ratio * 0.92 + header_bonus)
        reason = f"{matched} of {tested} evaluated values match {IDENTIFIER_LABELS[kind]}."
        scored.append(DetectionResult(kind, confidence, matched, tested, reason, multiple and kind == IdentifierType.UNIPROT))

    gene_headers = ("gene", "gene symbol", "gene names", "hgnc", "symbol")
    gene_matches, tested, _ = _ratio(values, _GENE_SYMBOL)
    gene_ratio = gene_matches / tested
    gene_hint = any(hint in header for hint in gene_headers)
    gene_confidence = min(0.82, gene_ratio * (0.55 if gene_hint else 0.28) + (0.22 if gene_hint else 0.0))
    scored.append(DetectionResult(
        IdentifierType.GENE_SYMBOL,
        gene_confidence,
        gene_matches,
        tested,
        f"{gene_matches} of {tested} values look like gene symbols; confirmation is recommended.",
    ))

    entrez_matches, tested, _ = _ratio(values, _INTEGER)
    entrez_ratio = entrez_matches / tested
    entrez_hint = any(hint in header for hint in ("entrez", "ncbi gene", "gene id"))
    entrez_confidence = entrez_ratio * (0.82 if entrez_hint else 0.18) + (0.12 if entrez_hint else 0.0)
    scored.append(DetectionResult(
        IdentifierType.ENTREZ,
        min(0.94, entrez_confidence),
        entrez_matches,
        tested,
        f"{entrez_matches} of {tested} values are integers; the header "
        + ("supports Entrez." if entrez_hint else "does not confirm they are Entrez IDs."),
    ))

    best = max(scored, key=lambda item: item.confidence)
    if best.confidence < 0.55:
        return DetectionResult(
            IdentifierType.UNKNOWN,
            best.confidence,
            best.matched_values,
            best.tested_values,
            "No identifier pattern reached sufficient confidence. " + best.reason,
            best.contains_multiple,
        )
    return best

