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

This release deliberately does not provide differential analysis or later biological modules.

## InterPro/Pfam annotation support

Database Manager retrieves compact release metadata exclusively from the official EMBL-EBI
InterPro API. Project workflows can then selectively acquire InterPro integrated entries and Pfam
signatures for an explicit collection of UniProt accessions. Responses are reused through a
release-aware global cache and materialized as immutable, hashed project annotation sets.

Raw responses, memberships, repeated and overlapping domain locations, entry types, and available
integration metadata are preserved. Once a set is Ready, all annotation loading and protein lookup
operations work offline.

The scientific backend analyzes unique UniProt accessions against an experimental target and
background without network access. It produces separate InterPro and Pfam frequency and enrichment
families, Benjamini-Hochberg correction, repeated-feature and overlap-preserving domain
architectures, plots, workbooks, immutable runs, and complete provenance.

The Analyses workspace provides the corresponding InterPro/Pfam scientific page. It can explicitly
build annotation sets, validate coverage, run analyses in the background, navigate proteins,
features and repeated locations, inspect architectures and historical runs, and export unchanged
tables, workbooks and plots. Analysis and historical loading remain offline after acquisition.

## MitoCarta3.0 database support

The Database Manager can install the official MitoCarta3.0 database for Homo sapiens as a reusable
local snapshot. It normalizes the human gene inventory, sub-mitochondrial compartment annotations,
and MitoPathways membership and hierarchy tables. After installation these database resources are
available offline.

The Analyses workspace provides the scientific MitoCarta3.0 workflow for human projects. It runs the
validated gene-level pipeline in the background and presents membership, overall enrichment,
sub-compartments, MitoPathways, bidirectional gene navigation, plots, unchanged artifact exports,
and immutable historical runs. Analysis remains fully offline; data installation is handled only by
Database Manager.

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

### STRING Homo sapiens local database

The Database Manager can acquire the official Homo sapiens STRING functional association network, physical network, protein metadata, and identifier aliases as an immutable local snapshot. Full evidence channels and raw scores are retained in an indexed SQLite database for reproducible offline mapping and queries. STRING data attribution and license: STRING Consortium, CC BY 4.0.

### Offline STRING network analysis

STRING analyses use the installed Homo sapiens SQLite snapshot to build experimental subnetworks from uniquely mapped UniProt seeds. Functional association and physical networks are analyzed in separate runs. PichAnalysis offers combined-score threshold presets of 150, 400, and 700, plus a custom 0–1000 value. Optional one-hop and two-hop expansion uses shortest-hop assignment; strict common direct neighbors must connect to every seed above the threshold. R igraph computes component membership, network degree, unweighted betweenness, component-local closeness, clustering, and explicit degree-based hubs on the run subnetwork. Full evidence channels, run-specific snapshot hashes, plots, and a workbook are preserved. Analysis and historical loading require no network access.

The Analyses → STRING page provides a local interactive network viewer for functional or physical runs. Researchers can select 150/400/700 confidence presets or a custom score, inspect Degree 1 and Degree 2 expansion, and use union or strict common direct-neighbor selection. Persisted summary, node metrics, hubs, components, support, and per-edge evidence remain available alongside plots. Viewer search, component and hop filters, node dragging, and image export change presentation only; they never recompute scientific results. Historical runs load offline from their own immutable run artifacts, including their original snapshot identity. Full node/edge CSVs, the existing workbook, and plot images can be exported without rerunning the analysis.

### Curated human Complex Portal database

The Database Manager downloads the official EMBL-EBI Complex Portal ComplexTab for Homo sapiens (NCBI taxonomy 9606) into immutable local snapshots. This database contains manually curated complexes only; predicted hu.MAP/MuSIC complexes are not included in the current PichAnalysis Complex Portal database. The original source file and its SHA-256 hash are retained alongside normalized direct participants, expanded protein membership, non-protein participants, nested-complex relationships, evidence, curator source, stoichiometry, and a searchable SQLite index. A source stoichiometry of 0 means unknown, not zero molecules. Alternative protein sets in the expanded source are preserved without being promoted to definite individual UniProt memberships. Local UniProt lookups can normalize valid isoform accessions to the canonical form for this module only, while retaining the input accession. Snapshot queries work offline. Complex coverage, presence, and enrichment are not inferred at this stage.

Complex Portal citation: Balu S. et al. (2025), *Complex Portal 2024: Predicted human complexes and enhanced visualisation tools for the comparison of orthologous and paralogous complexes*, Nucleic Acids Research, https://doi.org/10.1093/nar/gkae1085. Consult the Complex Portal terms at https://www.ebi.ac.uk/complexportal/about#license_privacy; an explicit data license is not asserted unless present in retrieved official metadata.
### Offline Complex Portal protein-component analysis

The Complex Portal analysis uses the installed manually curated human snapshot and unique canonical UniProt identifiers from the experimental protein catalog. Valid isoforms are normalized only for this module; original accessions and supporting source rows remain in the mapping outputs. The default background is all uniquely resolved canonical proteins in the experiment, including proteins with no curated complex membership. Presence/Absence-derived or manual targets can be selected without recomputing those upstream results.

Protein-component coverage is calculated over component groups from the official expanded representation. A single protein is one group; an alternative set such as `[P2,P3,P4]` is also one group and is covered when any permitted option is detected. Multiple detected alternatives still cover only one group, though each unique detected protein can count separately for frequency and enrichment. Stoichiometric copy counts do not multiply groups, and unknown stoichiometry remains unknown. Non-protein and nested participants are retained for audit but do not inflate the protein-component denominator.

The reported classes are complete protein-component coverage, partial protein-component coverage, no detected protein components, and not applicable for complexes without protein components. Complete protein-component coverage means every expected protein component group has at least one compatible protein identity in the target. It does **not** demonstrate physical assembly, simultaneous interaction, co-localization, correct stoichiometry, activity, or detection of RNA/chemical participants. Frequency counts unique target canonical proteins annotated as possible members. Enrichment uses one-sided Fisher/hypergeometric testing against the experimental background minus target, excludes complexes below the explicit minimum overlap, and applies Benjamini–Hochberg FDR over the eligible curated-complex family. Runs preserve their original snapshot, provenance, plots, workbook, and complete tables for offline historical inspection. No quantitative proteomics intensities are used to infer coverage.
The Analyses → Complexes page presents persisted protein-component coverage, alternative component groups, frequency and enrichment as distinct results. It keeps direct participants separate from expanded protein membership and shows non-protein participants, nested complexes, stoichiometry, evidence, and original isoform accessions in the complex detail and mapping views. Historical runs load offline from their own snapshot-bound artifacts, without recalculation. The current table, existing workbook, and existing PNG/PDF plots can be exported; the selected complex detail sections can also be exported as CSV.
### Managed Gene Ontology snapshot

The Database Manager can download immutable Homo sapiens Gene Ontology snapshots from the official GO `go-basic.obo` ontology and `HUMAN-uniprot.gaf.gz` annotation sources. Each snapshot retains the original downloads, source URLs, versions, SHA-256 hashes, normalized term/relation/annotation tables, and a local SQLite index. Mixed-taxon annotations are explicitly excluded and counted; `NOT` qualifiers remain auditable. Updates preserve the previous Ready snapshot on failure or cancellation. This managed data source supports mtDNA Evidence; it does not replace the existing scientific GO analysis based on R annotation packages.

### Offline mtDNA evidence index

The Database Manager offers a derived mtDNA Evidence index for Homo sapiens. It requires a Ready local MitoCarta snapshot and a Ready local Gene Ontology snapshot with normalized go_terms.csv, go_relations.csv, and go_annotations.csv tables; it does not download either dependency. The existing GO analysis uses local R annotation packages, which are not themselves a managed GO snapshot. The Database Manager can install this managed GO snapshot; the mtDNA Evidence card disables Build until both it and MitoCarta are Ready. The only newly retrieved source during a real build is the official NCBI RefSeq mitochondrial record NC_012920.1.

The index preserves separate MitoCarta pathway, mitochondrial-specific GO annotation, and mtDNA genomic-origin evidence records with source snapshot IDs, hashes, original terms, GO evidence codes and qualifiers. GO NOT annotations are audited but never positive evidence. MitoCarta membership alone, generic DNA annotations, and unrelated mtRNA metabolism are not mtDNA evidence. An mtDNA-encoded product is not automatically DNA-binding or nucleoid-associated. Source count is descriptive, not a global confidence score. Derived snapshots and their SQLite indexes remain queryable offline against their original source versions.
### Offline mtDNA Evidence analysis

Analyses → mtDNA Evidence uses only a Ready local derived mtDNA Evidence snapshot; it performs no network requests. The statistical unit is the unique gene-level evidence entity (`NCBI:<GeneID>` or a justified `SYMBOL:<SYMBOL>` fallback), never a proteomics row or individual evidence record. Original experimental rows and accessions remain in the mapping audit. Correctly mapped genes without mtDNA-related evidence remain in the default experimental background and are distinguished from unmapped or ambiguous rows.

R calculates descriptive category and source frequencies, independent-source agreement, and one-sided Fisher category enrichment against the experimental background minus target. Categories below the explicit overlap threshold are excluded; Benjamini–Hochberg adjustment applies to every tested category. `mtDNA_encoded` is a separate genomic-origin claim and does not imply nucleoid localization, DNA binding, replication or repair. Source count is descriptive, not a confidence score. Presence/Absence condition comparisons are descriptive only and do not establish recruitment or loss. Each run freezes the snapshot identity, source hashes, mapping, target/background, R scripts, outputs, plots and workbook. Historical runs load these persisted artifacts offline without consulting a newer active snapshot.
### Offline Proteomics QC scientific core

Proteomics QC uses only columns configured as Quantification and requires one compatible quantification family plus sample condition and replicate metadata. The feature unit is one original experimental row: protein groups are not exploded and duplicate biological identifiers are not collapsed. Python freezes the selected matrix, row identifiers, sample metadata, explicit transformation and zero policy, hashes, and R scripts in a run-specific directory. Historical runs load their own persisted outputs offline.

LFQ and raw intensities default to log2 of positive values with zeros treated as missing; no pseudocount is added. Spectral counts default to log2(x+1). Other quantitative values require an explicit zero policy. Detection is separate from transformation. Missingness, per-sample distributions, pairwise-complete Pearson and Spearman correlations, shared counts, detection-set Jaccard, and within-condition replicate summaries are reported without assigning a QC pass/fail score. Feature CV uses linear untransformed values, requires at least two finite detected measurements with positive mean, and never imputes. Centered, non-unit-scaled PCA uses complete-case features only; zero-variance features are excluded only from PCA. Unavailable metrics remain missing. Sample diagnostics and structural warnings never exclude samples automatically. Runs include CSV tables, PNG/PDF plots, and an Excel workbook; the scientific core does not perform normalization or network access.

### Proteomics QC interface and history

Analyses → Proteomics QC selects configured quantitative columns from one compatible family, shows condition and replicate metadata, and exposes the backend's transformation and zero-as-missing policies. A background worker runs the frozen Python/R QC pipeline while the interface remains responsive. The page displays persisted summary, sample metrics, detection, missingness, distributions, pairwise Pearson and Spearman correlations with shared counts, detection-set Jaccard, linear-scale CV, complete-case PCA, distances, diagnostics, structural warnings, and the R-generated plots. Metrics are descriptive: there is no imputation, QC score, automatic sample exclusion, normalization, or differential analysis. Users can filter tables, inspect samples, navigate graphs, export the current table or persisted workbook and PNG/PDF, and reopen historical runs offline. Historical results use their own frozen inputs and parameters even when the project's current input changes; incomplete runs are labeled instead of borrowing artifacts from another run.

### Offline differential matrix preparation (no differential model yet)

Differential preparation selects exactly two conditions and at least two chosen samples per condition from one compatible continuous quantitative family. Spectral counts require a separate count-based model; Other quantitative values require explicit confirmation that they are continuous. Each original experimental row remains one feature: protein groups are not exploded and duplicate biological identifiers are not collapsed. The run freezes the original selected matrix, sample and feature metadata, input hashes, condition order, and the future Condition A − Condition B contrast direction.

The pipeline records detection and missingness first, then classifies observed patterns and continuous-analysis eligibility before any transformation, normalization or imputation. The default requires at least two observed replicates in each condition and more than two total observations. Partially missing features can remain eligible; a feature undetected in a whole condition is excluded from the continuous matrix and preserved separately as a qualitative detection candidate. Sparse one-condition detections are distinguished from candidates meeting the selected threshold. These categories do not infer MAR, MNAR, biological absence or a detection-limit mechanism.

LFQ and raw intensities default to log2 of positive values with zero treated as missing and no pseudocount. Transformation may explicitly be disabled. Optional median centering uses observed transformed sample medians and the median of those medians as the reference; no missing values are filled for normalization. No imputation is the default and the prepared matrix may retain NA. Opt-in MinProb and QRILC are left-censored approaches; KNN is similarity-based. The application uses the official Bioconductor MsCoreUtils impute_matrix API and requires MsCoreUtils plus imputeLCMD for MinProb/QRILC or impute for KNN; packages are never installed during an analysis. Seed, method-specific parameters, imputed cells, unresolved cells and masks are preserved. Only partially missing eligible features are imputed; an entirely absent condition is never imputed and imputation cannot rescue an ineligible feature. Observed values must remain unchanged.

Each immutable offline run stores original, transformed, normalized and prepared CSV matrices; eligibility and qualitative audits; design preview; plots in PNG/PDF; an Excel workbook; R session information and provenance. Historical runs load their own frozen inputs. This stage does not fit limma models or calculate logFC, p-values or FDR and adds no differential-analysis GUI.

### Offline differential statistics with limma

Differential Statistics consumes one validated, frozen Differential Preparation run; it never reconstructs preparation from the current project table or re-imputes missing values. It copies the prepared matrix, sample and feature metadata, contrast, eligibility, qualitative candidates and imputation audit into an immutable child run, with parent IDs and hashes. The original experimental row remains the statistical feature, including protein groups and duplicate biological identifiers. A completed child can be reopened offline even if its parent is no longer available.

For two independent conditions, R fits a direct Intercept plus A_minus_B design using limma::lmFit and limma::eBayes; the coefficient is Condition A minus Condition B, with no contrasts.fit approximation or paired/block model. Missing prepared values remain missing, so residual degrees of freedom are feature-specific. Trend and robust empirical Bayes are explicit options and default to false. The limma version, design, variance moderation parameters, residual degrees of freedom, confidence intervals, moderated t, raw p-value, B statistic and AveExpr are retained. Robust moderation requires statmod, which is validated by the R environment.

A log2-positive preparation has log2 effects and may report log2FC, fold change and percent change. An untransformed preparation requires an explicit declaration of log2 or continuous scale; continuous effects are differences, never labelled fold changes. BH adjustment covers all successfully tested features before any effect-size classification. The FDR threshold and minimum absolute effect are separate; changing the latter never changes adjusted p-values. Qualitative detected-only candidates stay in a separate descriptive table without p-values or infinite fold changes. Imputation counts from the parent are attached per feature, but the fitted model does not propagate imputation uncertainty. Diagnostic volcano, MA, mean-comparison, p-value, FDR, effect, residual-df and top-effect plots are persisted as PNG/PDF; an imputation-fraction diagnostic is added only for an imputed parent. Plotting floors numerically zero adjusted p-values at machine epsilon without changing tabular values. The run also includes CSV tables, an Excel workbook and R provenance. No GUI, TREAT analysis or network access is added at this stage.

### Differential Analysis interface

Analyses → Differential Analysis offers a two-stage, offline workflow. Differential Preparation selects exactly two project conditions and explicit quantitative samples (at least two per condition), applies an explicit zero policy and log2-positive or no transformation, and determines feature eligibility before optional median centering and optional imputation. No imputation is the default; MinProb, QRILC, and KNN are opt-in and available only when their local R packages are installed. The application does not infer MAR/MNAR or exclude samples automatically. Detected-only features remain qualitative candidates outside the continuous model.

Differential Statistics consumes a Ready, hash-validated, frozen Preparation run. An untransformed parent requires an explicit log2 or continuous prepared-scale declaration. R limma fits the direct Condition A minus Condition B coefficient with lmFit and eBayes; trend and robust moderation are opt-in. BH-adjusted FDR is computed over tested features, while the effect threshold only classifies results. A continuous-scale effect is a difference, not fold change. The page displays persisted summary and full/tested/significant/untested results, qualitative candidates, model and imputation audit, and R-generated volcano/MA and diagnostic plots. Imputed values are visible per feature, but the model does not propagate imputation uncertainty.

Both stages have independent run histories. Historical tables, plots, provenance and workbooks are read from each run's frozen artifacts without contacting a newer project matrix or the network; completed Statistics runs can still be read when their parent Preparation directory is unavailable. Current tables, selected-feature details, the persisted workbook, and PNG/PDF plots can be exported. The Scripts / Logs page exposes the 12_differential_preparation.R and 13_differential_analysis.R entry points and per-run copies.
### Cross-module Explorer

Analyses → Cross-module Explorer searches a local, rebuildable index of persisted analysis results. It does not run a new scientific calculation or require network access. Search by feature ID, source row, original identifier, protein or gene identifier, or a module-specific target; select a source record to inspect its identity, run, snapshot and possible links. A selected result can be opened in its analysis page, and a cross-module summary can be exported as CSV.

An experimental feature is one original input row, not a protein or gene. Exact row navigation requires a compatible frozen input lineage; protein and gene identifiers provide biological bridges only. Protein groups are not exploded into quantitative observations, duplicate identifiers remain separate rows, ambiguous mappings are labelled, and UniProt isoforms are not globally canonicalized. The Explorer uses identifiers and outputs frozen in each run, preserves historical runs and database snapshot provenance, and never borrows another run's artifacts. Cross-lineage lookup is opt-in and shows biological matches only: a match from a different input lineage does not establish that two experimental features are the same.
Use **Explore across analyses...** from a selected result row in a loaded analysis. **Open in analysis** loads the explicitly selected persisted run, focuses an existing compatible control where safe, or displays the target identifier and run as a visible fallback. If several source records or target entities match, choose the intended row explicitly; the Explorer does not take the first match. Historical Mapping runs with only raw responses and no frozen mapping tables are reported as unavailable for safe resolution. The current/latest Mapping catalog is never substituted for an older run.
### Consolidated Reports

The Reports page builds an offline, presentation-only bundle from explicitly selected persisted runs and CSV/PNG artifacts. Select runs and artifacts, order sections, add optional user notes, and generate both HTML and PDF in the background. Reports preserve each source run ID, input lineage and historical database snapshot; the overview warns when selected runs have multiple input lineages. Tables show a bounded preview, with complete selected files copied as attachments when enabled. No analysis, remapping, statistical recalculation, network expansion or automated biological interpretation is performed.

Bundles are stored under `reports/consolidated/runs/<report_id>` with `report.html`, `report.pdf`, configuration, SHA-256 manifest, CSV source/artifact indexes, figures and optional attachments. History validates hashes before opening or exporting. Missing or altered artifacts are refused; incomplete staging directories remain marked Incomplete for audit. Historical Mapping runs without frozen mapping tables cannot borrow the latest catalog.
