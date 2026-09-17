# PichAnalysis

PichAnalysis is a cross-platform desktop foundation for reproducible proteomics
projects. It creates and opens projects, imports tabular data, previews columns,
helps interpret identifiers and quantitative columns, and checks a local R
installation. Scientific analyses are not implemented yet.

## Architecture

The PySide6 interface is separated from project and import logic. Python manages
the desktop workflow and data files; future scientific analyses will run through
a central `Rscript` adapter. Application-owned R scripts live in `r_scripts/`.

## Setup

Python 3.10 or newer is required. R is optional for the current data workflow,
but will be required by future scientific analyses.

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pichanalysis
```

macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pichanalysis
```

Run tests with `.\.venv\Scripts\python.exe -m pytest` on Windows or
`.venv/bin/python -m pytest` on macOS/Linux.

## Project layout

Each user-selected project contains `project.json`, original and processed input
folders, mapping and annotation folders, analyses, reports, logs, and per-run
scripts. Imported originals are copied without modification and internal data is
stored as UTF-8 CSV. A project can be moved and reopened because its internal
paths are relative to its root.

## Column configuration

After import, PichAnalysis suggests a role for each column, preliminary biological
identifier types, and condition/replicate labels for clearly named quantitative
columns. The user can correct every suggestion, choose one primary identifier,
and save the experimental design in `project.json` without changing the table.

Identifier detection is local and deterministic. It is based on column names and
patterns found in a bounded sample of values; it is **not validation against
UniProt, NCBI, Ensembl, HGNC, or any other external biological database**.

## Biological mapping and annotation

The project organism must be configured before mapping. Homo sapiens (9606), Mus
musculus (10090), and custom scientific-name/taxonomy-ID pairs are supported.
The configured primary identifier is mapped by R through the official UniProt ID
Mapping REST workflow. UniProt candidates are preserved for unique, ambiguous,
and unmapped inputs, checked against the selected organism, and supplemented with
NCBI Gene annotations through the current Datasets v2 `dataset_report` endpoints.
UniProt and NCBI fields remain separate; no generated or integrated summary is
created.

Local identifier suggestions from the import screen are pattern-based hints.
Online mapping is the distinct validation/annotation step and still requires the
user to review ambiguous results.

### R dependencies

R must provide `httr2`, `jsonlite`, `readr`, `DBI`, `RSQLite`, and `openxlsx`.
The application never installs packages automatically. `Rscript` is searched on
PATH, in a configured location, standard versioned Windows installations, and
common macOS locations.

### Cache, provenance, and outputs

A reusable SQLite cache is stored in the platform-appropriate user cache folder,
outside the source repository. Selecting “Atualizar anotações online” bypasses
existing entries and appends freshly retrieved data. Every project still receives
an exact snapshot in `mapping/raw/<run_id>/`, plus scripts, parameters, and R
session information under `scripts/runs/<run_id>_mapping_annotation/`.

Automatic outputs in `mapping/tables/` are:

- `id_mapping.csv`
- `protein_catalog.csv`
- `unmapped.csv`
- `ambiguous.csv`
- `protein_mapping.xlsx`, with separate Mapping, Protein catalog, Unmapped, and
  Ambiguous worksheets

## Presence/absence and reproducibility

The Presence / absence module uses only columns explicitly configured as
`Quantification`. Quantification types are never mixed automatically. Detection
is evaluated on raw imported values using a user-selected threshold; by default,
zero and missing values are absence and a value greater than zero is presence.
No normalization or imputation is performed.

Reproducibility can require either a minimum number or a minimum fraction of
replicates. For two selected conditions, the strict classes include each
condition's `specific`, `Shared`, and `Sporadic`. Specific requires reproducible
detection in one condition and no detection in any replicate of the other.
`Predominantly <condition>` is a separate exploratory label and is never treated
as specific. Three or more conditions use an explicit reproducible-detection
pattern instead of invented pairwise classes.

Each run preserves tables, metadata, parameters, R session information, and the
scripts used. Latest outputs include the binary presence matrix, original
quantitative values, condition-level detection, classification, condition
summary, and an Excel workbook. R also generates protein counts per replicate,
an UpSet-style intersection plot, a bounded binary heatmap, and—for two
conditions—a detection-fraction scatter plot. This module does not perform
differential statistical testing.

This release deliberately does not provide identifier mapping, statistics, GO,
KEGG, STRING, differential analysis, or other biological analyses.
