"""Canonical application release version."""

__version__ = "0.1.0"


def windows_file_version() -> str:
    return f"{__version__}.0"
