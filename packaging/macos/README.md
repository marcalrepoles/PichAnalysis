# macOS release preparation

Build the `.app` on a real macOS host with an explicitly selected architecture
(`arm64` or `x86_64`). This Windows task does not produce or validate a macOS
binary, universal2 bundle or DMG. Use a clean Python 3.11+ virtual environment,
install the project plus `packaging/requirements-build.txt`, run Python and R
gates, then execute `build_macos.sh arm64` or `build_macos.sh x86_64`.

The macOS bundle includes application-owned `r_scripts/` using the existing
frozen resource resolver. R remains external, with discovery through the
standard R framework and Homebrew locations. Scientific databases are not
bundled. Validate startup, report generation, R execution, installation and
uninstallation on macOS before calling it a release.

Developer ID signing, notarization and a DMG are future release steps. Do not
claim they occurred without actual Apple credentials and macOS validation.
