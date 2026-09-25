from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

from .column_mapping import suggest_columns
from .project import Project, ProjectError


class ImportError(Exception):
    """Raised when tabular input cannot be imported safely."""


SUPPORTED_SUFFIXES = {".xlsx", ".csv", ".tsv"}


@dataclass
class ImportResult:
    dataframe: pd.DataFrame
    original_path: Path
    processed_path: Path
    source_format: str
    sheet: str | None

    @property
    def column_metadata(self) -> list[dict[str, object]]:
        return [
            {
                "name": str(column),
                "dtype": str(self.dataframe[column].dtype),
                "non_null": int(self.dataframe[column].notna().sum()),
                "missing": int(self.dataframe[column].isna().sum()),
                "role": None,
                "identifier_type": None,
                "condition": None,
                "replicate": None,
                "quantification_type": None,
            }
            for column in self.dataframe.columns
        ]


def xlsx_sheets(path: Path) -> list[str]:
    try:
        with pd.ExcelFile(path, engine="openpyxl") as workbook:
            return list(workbook.sheet_names)
    except (OSError, ValueError, ImportError, Exception) as error:
        raise ImportError(f"Could not read the XLSX file: {error}") from error


def _detect_delimiter(path: Path, suffix: str) -> str:
    expected = "\t" if suffix == ".tsv" else ","
    try:
        sample = path.read_text(encoding="utf-8-sig", errors="strict")[:65536]
    except (OSError, UnicodeError) as error:
        raise ImportError(f"Could not read the UTF-8 text file: {error}") from error
    if not sample.strip():
        raise ImportError("The file is empty.")
    try:
        detected = csv.Sniffer().sniff(sample, delimiters=",;\t").delimiter
    except csv.Error:
        if expected not in sample.partition("\n")[0]:
            raise ImportError("Could not determine the file delimiter.")
        detected = expected
    first_line = sample.partition("\n")[0]
    candidates = [separator for separator in (",", ";", "\t") if separator in first_line]
    if len(candidates) > 1 and detected != expected:
        raise ImportError("The file delimiter is ambiguous; review the CSV/TSV file before importing.")
    return detected


def _validate_header(names) -> None:
    normalized = [str(name).strip() if name is not None else "" for name in names]
    if not normalized or any(not name for name in normalized):
        raise ImportError("Column names must not be empty.")
    if len(set(normalized)) != len(normalized):
        raise ImportError("Duplicate column names are not allowed.")


def _validate_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or frame.dropna(how="all").empty:
        raise ImportError("The input has no usable data rows.")
    return frame

def read_table(path: Path, sheet: str | None = None) -> tuple[pd.DataFrame, str | None]:
    path = Path(path)
    if not path.is_file():
        raise ImportError("The selected file does not exist.")
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_SUFFIXES:
        raise ImportError("Unsupported format. Use XLSX, CSV, or TSV.")
    try:
        if suffix == ".xlsx":
            sheets = xlsx_sheets(path)
            selected = sheet or (sheets[0] if len(sheets) == 1 else None)
            if selected is None:
                raise ImportError("Select a worksheet to continue.")
            if selected not in sheets:
                raise ImportError("The selected worksheet does not exist.")
            workbook = load_workbook(path, read_only=True, data_only=True)
            try:
                _validate_header(next(workbook[selected].values, ()))
            finally:
                workbook.close()
            return _validate_rows(pd.read_excel(path, sheet_name=selected, engine="openpyxl")), selected
        delimiter = _detect_delimiter(path, suffix)
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            _validate_header(next(csv.reader(stream, delimiter=delimiter), ()))
        return _validate_rows(pd.read_csv(path, sep=delimiter, encoding="utf-8-sig")), None
    except ImportError:
        raise
    except Exception as error:
        raise ImportError(f"Could not import the table: {error}") from error


def import_into_project(project: Project, source: Path, sheet: str | None = None) -> ImportResult:
    source = Path(source)
    dataframe, selected_sheet = read_table(source, sheet)
    try:
        original = project.copy_original(source)
        processed = project.root / "input" / "processed" / f"{original.stem}.csv"
        counter = 2
        while processed.exists():
            processed = processed.with_name(f"{original.stem}_{counter}.csv")
            counter += 1
        dataframe.to_csv(processed, index=False, encoding="utf-8")
    except (OSError, ProjectError) as error:
        raise ImportError(str(error)) from error
    result = ImportResult(dataframe, original, processed, source.suffix.lower()[1:], selected_sheet)
    project.config["input"] = {
        "original_file": project.relative(original),
        "processed_file": project.relative(processed),
        "format": result.source_format,
        "sheet": selected_sheet,
        "rows": int(len(dataframe.index)),
        "columns": int(len(dataframe.columns)),
        "column_metadata": result.column_metadata,
    }
    project.config["columns"] = suggest_columns(dataframe)
    project.config["column_configuration"] = {
        "status": "suggested",
        "errors": ["Review and save the suggested configuration."],
        "warnings": [],
    }
    project.save()
    return result
