# Runtime dependencies

This inventory reflects the source tree at the Prompt 18 hardening review. The
application never installs Python or R packages during an analysis.

## Python

`pyproject.toml` declares PySide6, pandas, openpyxl, xlrd, and platformdirs for
runtime. Pytest is development-only. Standard-library `urllib` provides the
database download clients; no additional Python HTTP package is required.

## R

An independently installed `Rscript` is required for scientific scripts, not
for opening the application or browsing persisted results. The release
preflight reports package name, required-by module, installation state, and
version without installing anything.

| Package | Required by |
| --- | --- |
| jsonlite, readr, openxlsx | Mapping and multiple scientific modules |
| httr2, DBI, RSQLite | Identifier mapping and cache |
| ggplot2 | Scientific plots |
| AnnotationDbi, GO.db | GO core |
| org.Hs.eg.db, org.Mm.eg.db | GO, conditional on the project organism |
| xml2 | KEGG KGML |
| igraph | STRING |
| limma | Differential Statistics |
| MsCoreUtils | Opt-in Differential Preparation imputation |
| imputeLCMD | Opt-in MinProb/QRILC |
| impute | Opt-in KNN |
| statmod | Opt-in robust limma moderation |

R base/recommended packages (`stats`, `utils`, `graphics`, `grDevices`) are not
separate installer downloads. The local R test suites have additional test
requirements but are not application runtime dependencies.

## Network and local data

Identifier mapping uses NCBI/UniProt-backed R helpers when configured to query
online. Database acquisition uses official GO, KEGG, Reactome, MitoCarta,
STRING, Complex Portal, and NCBI sources. Selective InterPro retrieval uses
the InterPro API. Scientific analyses with already installed local snapshots,
Cross-module Explorer, and Consolidated Reports use persisted local data. The
installer must not silently download scientific databases or R packages.

Database licensing and redistribution terms must be verified before bundling
third-party datasets; this document does not grant redistribution rights.
