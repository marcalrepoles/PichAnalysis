# Final hardening review

Prompt 18 is a defect-focused engineering audit. Scientific methods are frozen.
The issue table records reproduced behavior rather than speculative cleanup.

| Area | Severity | Reproduction | Correction | Regression | Status |
| --- | --- | --- | --- | --- | --- |
| R discovery | Major | An explicit Rscript path lost to a different executable on `PATH` | Honor a configured path first; an invalid explicit path reports missing | `test_explicit_rscript_configuration_wins_over_path` | Fixed |
| CSV/TSV/XLSX import | Major | Duplicate headers were silently renamed by pandas | Validate original headers before pandas parsing | `test_duplicate_input_columns_are_rejected`, `test_duplicate_xlsx_columns_are_rejected` | Fixed |
| Empty import | Minor | Header-only and all-blank datasets were accepted | Reject data with no usable rows | `test_header_only_or_blank_rows_are_rejected` | Fixed |
| Project input reload | Major | Missing processed/original input was logged but not shown to the user | Validate both paths and show missing-artifact message in Data page | `smoke_application_hardening.py` | Fixed |
| App shutdown | Major | On Windows, the project log remained open after closing MainWindow, preventing project relocation/removal | Close and detach log handlers on accepted close | `smoke_application_hardening.py` | Fixed |
| R script resources | Blocker for packaging | Scientific modules assumed `parents[3]/r_scripts`, valid only in source layout | Resolve source/frozen/explicit resource roots; validate missing scripts | `test_resource_override_and_different_cwd`, application smoke | Fixed for the audited layouts; installer must include scripts |
| Consolidated Report PDF | Blocker | Selected PNGs rendered as broken-image icons in PDF despite intact HTML | Register local image resources with QTextDocument before printing | `test_explicit_offline_report_and_history` renders PDF and checks original-blue figure pixels | Fixed |
| Offline smoke imports | Minor | Four smoke scripts imported an unrelated `tests` package and stopped before exercising code | Use sibling smoke imports without changing pytest package layout | Four smoke scripts rerun in new processes | Fixed |

## Coverage and boundaries

Project creation/opening, Unicode/space paths, relocation, missing/corrupt or
unsupported metadata, CSV/TSV/XLSX import, explicit sheet selection, R
discovery, main navigation, empty histories, missing-R startup, clean shutdown,
and CWD independence were covered by unit tests and the fresh-process smoke.
The existing Python and R suites cover database snapshot handling, scientific
run loaders, historical isolation, workers, Cross-module Explorer, and reports.
No new scientific model, interpretation, database integration, or installer was
introduced. See `runtime-dependencies.md` for the audited package/network
inventory and `release-resources.md` for future installer inputs.

The macOS review is static: `pathlib` and Qt file URLs provide platform-neutral
paths; R discovery includes standard framework/Homebrew locations. Runtime
execution on macOS remains to be verified during Prompt 19. Third-party
database redistribution terms also remain to be verified before bundling.

## Prompt 18-Fix language review

Inspected `src/pichanalysis/ui/`, `src/pichanalysis/core/`, and user-visible
messages/plot labels in `r_scripts/`. Portuguese text found in project/import
errors, column and identifier detection, Mapping/GO/Presence readiness,
Data/Analyses/GO/Presence/Project/Scripts pages, session logs, and R errors
and plot labels. These displayed strings were translated to English without
changing identifiers, parameters, defaults, thresholds, or scientific methods.
The expanded `test_ui_english.py` checks concrete reviewed phrases; a manual
re-search for accented text and the identified ASCII-only phrases returned no
runtime matches. Portuguese intentionally remains in historical documentation,
test fixture names, and tests of Unicode paths; user-supplied data and names
are not translated.

## Remaining known issues

- Blocker: None known
- Major: None known
- Minor:
  - macOS execution not yet validated on macOS hardware
  - full real STRING acquisition smoke not executed during final hardening

## Final validation

- Python: 280 passed, 0 failed, 0 skipped after Prompt 18-Fix.
- R: 13/13 scientific suites passed again after text-only R translations.
- Offline smokes: 12/12 passed in separate processes.
- Real-data acquisition smokes: InterPro, KEGG, MitoCarta, Complex Portal, Reactome, and mtDNA Evidence passed in temporary directories. The mtDNA smoke exercised real GO and MitoCarta snapshots.
- Full STRING real-data smoke: not run; its offline fixtures, Python regressions, and R scientific suite passed.
- `python -m compileall -q src`: passed.
- Report PDF: first and figure/table pages rendered and visually reviewed; QtPdf regression checks original figure pixels.
