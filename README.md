# PichAnalysis

PichAnalysis is a cross-platform desktop application for reproducible proteomics
projects. It creates and opens projects, imports tabular data, maps identifiers,
and provides local Presence/Absence, Gene Ontology, KEGG, and Reactome analyses.

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

## Gene Ontology

The GO module keeps individual protein-to-GO annotation, descriptive term
frequency, and statistical enrichment as separate products. Frequency answers
how often a term occurs; enrichment compares the target with an explicit
experimental background using a hypergeometric test. The default background is
mapped proteins from the experiment, never the whole genome implicitly.

BP, MF, and CC are processed separately. Direct OrgDb annotations retain
evidence codes and may include all evidence, exclude IEA, or retain only a
centralized set of experimental codes. Enrichment reports raw p-values,
Benjamini-Hochberg FDR, GeneRatio, BgRatio, counts, and associated genes.
Optional Jaccard redundancy simplification never replaces complete results.

Local support uses `AnnotationDbi`, `GO.db`, `org.Hs.eg.db`, and `org.Mm.eg.db`
for human and mouse. Runs record R, GO.db, and OrgDb versions and preserve long
annotations, unannotated proteins, per-ontology tables, workbook, plots,
metadata, and immutable script snapshots under `analyses/GO/`.

This release deliberately does not provide MitoCarta, STRING, differential
analysis, or other later biological modules.

## Reactome analysis

The Reactome analysis page uses an active local Homo sapiens Core Data snapshot
installed through Database Manager. Analysis is fully offline: UniProt-first and
NCBI Gene fallback mapping, pathway membership, hierarchy, frequency, and
hypergeometric enrichment with Benjamini-Hochberg FDR are computed by the
versioned R pipeline. The experimental Reactome-mapped set is the default
background.

The GUI provides target and background configuration, advanced statistical
options, local pathway-to-protein and protein-to-pathway navigation, generated
plots, exports, and immutable historical run loading. Ambiguous and unmapped
entities remain visible. The page does not download Reactome data or provide a
local diagram viewer. The optional official `diagrams.png.tgz` component is installed as a new atomic
snapshot, independently of Reactome Core Data. The Database Manager shows its status and supports
download, update, progress, and cancellation. Pathway results can open the original local PNG with
fit/actual-size/zoom controls and export an unchanged copy. Reactome diagrams are attributed to
Reactome and licensed under CC BY 4.0; the tabular Core Data remains CC0.

## Database Manager

The **Database Manager** page is available without an open analysis project. It installs or updates a local Homo sapiens (`hsa`, NCBI Taxonomy ID 9606) KEGG snapshot containing normalized core tables and, optionally, pathway entries, KGML files, and PNG images.

Data is stored in the operating system's per-user application-data directory, outside projects and the source tree. Each update uses a new snapshot. `active.json` changes atomically only after completion, so a failed or cancelled update cannot replace a working database. Incomplete downloads can resume and skip valid completed files.

KEGG REST access is centrally limited to at most three calls per second with bounded retries. The application displays an academic-use notice before the first download. Review <https://www.kegg.jp/kegg/rest/> and <https://www.kegg.jp/kegg/legal.html>.

For a small online smoke check, use a temporary root and call `DatabaseManager.download(..., pathway_subset=["hsa00010", "hsa00020"])`. Automated tests use fake responses and never write to the production database directory.

## KEGG Pathways analysis

KEGG Pathways requires an active, pathway-analysis-compatible Homo sapiens snapshot installed separately through Database Manager. KEGG data are not bundled with PichAnalysis. New snapshots include official local NCBI GeneID and UniProt conversion tables; historical snapshots remain untouched and may be structurally valid but incompatible with this analysis.

The analysis runs entirely offline in R. Stable NCBI Gene IDs are mapped first, followed by UniProt accessions. Users can exclude ambiguous mappings (the default) or include every mapped candidate for explicitly exploratory set analysis. Statistical units are deduplicated KEGG genes, while all originating experimental entities remain traceable.

Pathway frequency is descriptive and reports how many selected genes occur in each pathway. Pathway enrichment is a separate hypergeometric analysis against an explicit experimental background, with Benjamini-Hochberg FDR correction. The default background is the KEGG-mapped portion of the experiment, not the complete human genome.

The local Pathway Viewer can display installed PNG maps. KGML is required for highlighting target genes because the R analysis uses KGML membership and coordinates to identify hit nodes; Python only renders those calculated coordinates. Missing images or KGML do not prevent frequency and enrichment analysis.
