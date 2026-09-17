from __future__ import annotations

from dataclasses import dataclass
from .project import Project


@dataclass(frozen=True)
class Organism:
    name: str
    tax_id: str


COMMON_ORGANISMS = (Organism("Homo sapiens", "9606"), Organism("Mus musculus", "10090"))


def get_organism(project: Project) -> Organism | None:
    name = str(project.config.get("organism_name") or "").strip()
    tax_id = str(project.config.get("organism_tax_id") or "").strip()
    return Organism(name, tax_id) if name and tax_id else None


def has_biological_results(project: Project) -> bool:
    return any((project.root / "mapping" / "tables").glob("protein_catalog*"))


def set_organism(project: Project, name: str, tax_id: str) -> None:
    clean_name, clean_tax_id = name.strip(), tax_id.strip()
    if not clean_name:
        raise ValueError("Enter the organism's scientific name.")
    if not clean_tax_id.isdigit() or int(clean_tax_id) <= 0:
        raise ValueError("Enter a valid numeric NCBI Taxonomy ID.")
    project.config["organism_name"] = clean_name
    project.config["organism_tax_id"] = clean_tax_id
    project.save()
