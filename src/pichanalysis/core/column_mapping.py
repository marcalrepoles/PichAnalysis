from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any

import pandas as pd

from .identifier_detection import DetectionResult, IdentifierType, detect_identifier
from .project import Project


class ColumnRole(str, Enum):
    IGNORE = "ignore"
    IDENTIFIER = "identifier"
    QUANTIFICATION = "quantification"
    ANNOTATION = "annotation"
    OTHER = "other"


class QuantificationType(str, Enum):
    RAW_INTENSITY = "raw_intensity"
    LFQ_INTENSITY = "lfq_intensity"
    SPECTRAL_COUNT = "spectral_count"
    OTHER = "other_quantitative"
    UNKNOWN = "unknown"


ROLE_LABELS = {
    ColumnRole.IGNORE: "Ignore",
    ColumnRole.IDENTIFIER: "Identifier",
    ColumnRole.QUANTIFICATION: "Quantification",
    ColumnRole.ANNOTATION: "Annotation",
    ColumnRole.OTHER: "Other",
}
QUANTIFICATION_LABELS = {
    QuantificationType.RAW_INTENSITY: "Raw intensity",
    QuantificationType.LFQ_INTENSITY: "LFQ intensity",
    QuantificationType.SPECTRAL_COUNT: "Spectral count",
    QuantificationType.OTHER: "Other quantitative",
    QuantificationType.UNKNOWN: "Unknown",
}


@dataclass
class ColumnConfig:
    role: str = ColumnRole.OTHER.value
    identifier_type: str = IdentifierType.UNKNOWN.value
    primary_identifier: bool = False
    condition: str = ""
    replicate: str = ""
    quantification_type: str = QuantificationType.UNKNOWN.value
    detection: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MappingValidation:
    errors: tuple[str, ...]
    warnings: tuple[str, ...] = ()

    @property
    def valid(self) -> bool:
        return not self.errors


def _quantification_type(name: str) -> QuantificationType:
    lowered = name.casefold()
    if "lfq" in lowered:
        return QuantificationType.LFQ_INTENSITY
    if "spectral count" in lowered or "spectral_count" in lowered:
        return QuantificationType.SPECTRAL_COUNT
    if "intensity" in lowered:
        return QuantificationType.RAW_INTENSITY
    return QuantificationType.UNKNOWN


def _condition_replicate(name: str) -> tuple[str, str]:
    cleaned = re.sub(
        r"(?i)\b(?:lfq\s+intensity|raw\s+intensity|intensity|abundance|quantity|spectral\s+count)\b",
        "",
        name,
    ).strip(" _-.")
    match = re.fullmatch(r"(.+?)[ _.-]+(?:rep(?:licate)?[ _.-]*)?(\d+)", cleaned, re.IGNORECASE)
    if not match:
        return "", ""
    condition = match.group(1).strip(" _-.")
    return (condition, match.group(2)) if condition else ("", "")


def suggest_column(name: str, series: pd.Series) -> ColumnConfig:
    detection = detect_identifier(name, series)
    detection_data = {
        "detected_type": detection.detected_type.value,
        "confidence": round(detection.confidence, 4),
        "matched_values": detection.matched_values,
        "tested_values": detection.tested_values,
        "reason": detection.reason,
        "contains_multiple": detection.contains_multiple,
    }
    lowered = name.casefold()
    quantitative_hint = any(token in lowered for token in ("intensity", "lfq", "abundance", "quantity", "spectral count"))
    annotation_hint = any(token in lowered for token in ("description", "function", "localization", "coverage", "mass", "molecular weight", "go term"))

    if detection.detected_type != IdentifierType.UNKNOWN and detection.confidence >= 0.55:
        return ColumnConfig(role=ColumnRole.IDENTIFIER.value, identifier_type=detection.detected_type.value, detection=detection_data)
    if quantitative_hint and pd.api.types.is_numeric_dtype(series):
        condition, replicate = _condition_replicate(name)
        return ColumnConfig(
            role=ColumnRole.QUANTIFICATION.value,
            condition=condition,
            replicate=replicate,
            quantification_type=_quantification_type(name).value,
            detection=detection_data,
        )
    if annotation_hint:
        return ColumnConfig(role=ColumnRole.ANNOTATION.value, detection=detection_data)
    return ColumnConfig(role=ColumnRole.OTHER.value, detection=detection_data)


def suggest_columns(frame: pd.DataFrame) -> dict[str, dict[str, Any]]:
    suggestions = {str(name): suggest_column(str(name), frame[name]).to_dict() for name in frame.columns}
    identifiers = [name for name, config in suggestions.items() if config["role"] == ColumnRole.IDENTIFIER.value]
    if identifiers:
        best = max(identifiers, key=lambda name: suggestions[name]["detection"]["confidence"])
        suggestions[best]["primary_identifier"] = True
    return suggestions


def validate_mapping(columns: dict[str, dict[str, Any]]) -> MappingValidation:
    errors: list[str] = []
    warnings: list[str] = []
    identifiers = [(name, item) for name, item in columns.items() if item.get("role") == ColumnRole.IDENTIFIER.value]
    primary = [name for name, item in identifiers if item.get("primary_identifier")]
    if not primary:
        errors.append("No primary identifier selected.")
    elif len(primary) > 1:
        errors.append("Select only one primary identifier.")
    for name, item in identifiers:
        if item.get("identifier_type", IdentifierType.UNKNOWN.value) == IdentifierType.UNKNOWN.value:
            errors.append(f"Column {name} is marked as Identifier but has no defined type.")
        if item.get("detection", {}).get("contains_multiple"):
            warnings.append(f"Column {name} contains cells with multiple identifiers.")
    combinations: dict[tuple[str, str], str] = {}
    for name, item in columns.items():
        if item.get("role") != ColumnRole.QUANTIFICATION.value:
            continue
        condition = str(item.get("condition", "")).strip()
        replicate = str(item.get("replicate", "")).strip()
        if not condition:
            errors.append(f"Column {name} is Quantification but has no condition.")
        if condition and replicate:
            key = (condition.casefold(), replicate.casefold())
            if key in combinations:
                errors.append(
                    f"Columns {combinations[key]} and {name} use the same condition/replicate combination."
                )
            else:
                combinations[key] = name
    return MappingValidation(tuple(errors), tuple(warnings))


def save_mapping(project: Project, columns: dict[str, dict[str, Any]]) -> MappingValidation:
    validation = validate_mapping(columns)
    project.config["columns"] = columns
    project.config["column_configuration"] = {
        "status": "valid" if validation.valid else "invalid",
        "errors": list(validation.errors),
        "warnings": list(validation.warnings),
    }
    project.save()
    return validation


def load_mapping(project: Project) -> dict[str, dict[str, Any]]:
    value = project.config.get("columns", {})
    return value if isinstance(value, dict) else {}


def experimental_design(project: Project) -> dict[str, Any]:
    columns = load_mapping(project)
    primary_name = next((name for name, item in columns.items() if item.get("primary_identifier")), None)
    primary = columns.get(primary_name, {}) if primary_name else {}
    quantitative = [
        {
            "column": name,
            "condition": item.get("condition", ""),
            "replicate": item.get("replicate", ""),
            "quantification_type": item.get("quantification_type", QuantificationType.UNKNOWN.value),
        }
        for name, item in columns.items()
        if item.get("role") == ColumnRole.QUANTIFICATION.value
    ]
    return {
        "primary_identifier_column": primary_name,
        "primary_identifier_type": primary.get("identifier_type"),
        "quantification_columns": quantitative,
        "conditions": sorted({item["condition"] for item in quantitative if item["condition"]}),
        "replicates": sorted({item["replicate"] for item in quantitative if item["replicate"]}),
    }
