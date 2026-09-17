# PichAnalysis

PichAnalysis is a cross-platform desktop foundation for reproducible proteomics
projects. This first version creates and opens projects, imports tabular data,
previews columns, and checks a local R installation. Scientific analyses are not
implemented yet.

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

This release deliberately does not provide identifier mapping, statistics, GO,
KEGG, STRING, differential analysis, or other biological analyses.

