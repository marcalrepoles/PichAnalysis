# Install PichAnalysis on Windows

PichAnalysis 0.1.0 is distributed as a Windows x86_64 per-user installer.
Run `PichAnalysis-0.1.0-Windows-x86_64-Setup.exe`. The default location is
`%LOCALAPPDATA%\Programs\PichAnalysis`; administrator rights are not
normally required. Launch it from the Start Menu. Python and PySide6 are
included in the installer; do not install Python separately for this release.

R is external. The release was validated with R 4.6.0, but no minimum R
version is asserted. Install a compatible R separately before running
R-dependent analyses. PichAnalysis starts and permits non-R navigation even
when R is absent. The read-only release preflight checks Rscript, required R
packages, application R scripts, Python imports and writable paths. See
`runtime-dependencies.md` for the package inventory. The installer does not
install R or R packages, and never changes an existing R library. Install
missing packages manually using your normal R/Bioconductor procedure.

Scientific databases are not included. Open Database Manager to acquire
the datasets needed for your selected module. Database acquisition and
identifier-mapping network calls occur only through the existing user-driven
workflows; the installer itself performs no scientific downloads.

Projects and user-selected database folders are separate from the install
directory. Back up projects independently. Uninstall PichAnalysis through
Windows Installed Apps or its Start Menu uninstall entry; this removes the
application, not projects or user data outside the installation directory.
The Windows build is unsigned, so Windows may show a publisher warning.
