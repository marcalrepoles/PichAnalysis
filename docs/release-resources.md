# Resources required by a future installer

Prompt 18 does not build an installer. Prompt 19 must include:

- The `pichanalysis` Python package and declared Python runtime dependencies.
- `r_scripts/` with entry scripts `00_runtime_test.R` through
  `13_differential_analysis.R` and their `lib/` helpers, preserving relative
  paths for R `source(...)` calls.
- Existing application icons/assets, if added to a future package; none are
  currently required by the Python UI startup path.
- Project templates/resources only if later introduced; project structure is
  currently created by Python code.
- License/notice material for bundled code and any third-party data, after
  distribution rights are verified. Local scientific databases are not assumed
  to be bundled.

At runtime, scripts are resolved from `PICHANALYSIS_RESOURCE_ROOT/r_scripts`,
the frozen bundle root (`sys._MEIPASS/r_scripts`), beside a frozen executable,
or the source-tree `r_scripts/`. The process current working directory is not
used. `run_preflight()` verifies the critical scripts, Rscript, packages,
Python imports, and writable temp/project directories. Missing resources are
reported; no package installation or database download occurs.
