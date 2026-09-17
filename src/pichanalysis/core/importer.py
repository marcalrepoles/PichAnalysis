from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

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
        raise ImportError(f"Não foi possível ler o arquivo XLSX: {error}") from error


def _detect_delimiter(path: Path, suffix: str) -> str:
    expected = "\t" if suffix == ".tsv" else ","
    try:
        sample = path.read_text(encoding="utf-8-sig", errors="strict")[:65536]
    except (OSError, UnicodeError) as error:
        raise ImportError(f"Não foi possível ler o arquivo de texto em UTF-8: {error}") from error
    if not sample.strip():
        raise ImportError("O arquivo está vazio.")
    try:
        detected = csv.Sniffer().sniff(sample, delimiters=",;\t").delimiter
    except csv.Error:
        if expected not in sample.partition("\n")[0]:
            raise ImportError("Não foi possível determinar o separador do arquivo.")
        detected = expected
    first_line = sample.partition("\n")[0]
    candidates = [separator for separator in (",", ";", "\t") if separator in first_line]
    if len(candidates) > 1 and detected != expected:
        raise ImportError("O separador do arquivo é ambíguo; revise CSV/TSV antes de importar.")
    return detected


def read_table(path: Path, sheet: str | None = None) -> tuple[pd.DataFrame, str | None]:
    path = Path(path)
    if not path.is_file():
        raise ImportError("O arquivo selecionado não existe.")
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_SUFFIXES:
        raise ImportError("Formato não suportado. Use XLSX, CSV ou TSV.")
    try:
        if suffix == ".xlsx":
            sheets = xlsx_sheets(path)
            selected = sheet or (sheets[0] if len(sheets) == 1 else None)
            if selected is None:
                raise ImportError("Selecione uma planilha para continuar.")
            if selected not in sheets:
                raise ImportError("A planilha selecionada não existe.")
            return pd.read_excel(path, sheet_name=selected, engine="openpyxl"), selected
        delimiter = _detect_delimiter(path, suffix)
        return pd.read_csv(path, sep=delimiter, encoding="utf-8-sig"), None
    except ImportError:
        raise
    except Exception as error:
        raise ImportError(f"Não foi possível importar a tabela: {error}") from error


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
    project.save()
    return result

