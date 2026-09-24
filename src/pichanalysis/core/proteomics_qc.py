"""Offline, reproducible quantitative proteomics QC run preparation and loading.

Python validates and freezes inputs. Scientific metrics and plots are computed in R.
"""
from __future__ import annotations

import hashlib
import json
import math
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .column_mapping import ColumnRole, QuantificationType, experimental_design
from .r_runtime import RRuntime


class ProteomicsQCError(RuntimeError):
    code = "proteomics_qc_error"
    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


class NoQuantitativeColumnsError(ProteomicsQCError): code = "no_quantitative_columns"
class TooFewSamplesError(ProteomicsQCError): code = "fewer_than_two_samples"
class MixedQuantificationTypesError(ProteomicsQCError): code = "mixed_quantification_types"
class MissingConditionError(ProteomicsQCError): code = "missing_condition_metadata"
class MissingReplicateError(ProteomicsQCError): code = "missing_replicate_metadata"
class InvalidNumericValueError(ProteomicsQCError): code = "invalid_numeric_value"
class NegativeValueError(ProteomicsQCError): code = "negative_value_incompatible"
class UnsupportedTransformationError(ProteomicsQCError): code = "unsupported_transformation"
class NoUsableQuantitativeDataError(ProteomicsQCError): code = "no_usable_quantitative_data"
class QCRExecutionError(ProteomicsQCError): code = "r_execution_failure"
class MissingQCOutputError(ProteomicsQCError): code = "missing_required_output"


DEFAULT_TRANSFORMATION = {
    QuantificationType.LFQ_INTENSITY.value: "log2_positive",
    QuantificationType.RAW_INTENSITY.value: "log2_positive",
    QuantificationType.SPECTRAL_COUNT.value: "log2p1",
    QuantificationType.OTHER.value: "none",
}
DEFAULT_ZERO_MISSING = {
    QuantificationType.LFQ_INTENSITY.value: True,
    QuantificationType.RAW_INTENSITY.value: True,
    QuantificationType.SPECTRAL_COUNT.value: True,
}
TRANSFORMATIONS = {"log2_positive", "log2p1", "none"}


@dataclass(frozen=True)
class QCParameters:
    selected_columns: tuple[str, ...]
    transformation: str | None = None
    zero_is_missing: bool | None = None


TABLES = {
    "sample_detection": "detection/sample_detection.csv",
    "feature_detection": "detection/feature_detection.csv",
    "feature_condition_detection": "detection/feature_condition_detection.csv",
    "replicate_consistency": "detection/replicate_consistency.csv",
    "sample_missingness": "missingness/sample_missingness.csv",
    "feature_missingness": "missingness/feature_missingness.csv",
    "condition_missingness": "missingness/condition_missingness.csv",
    "sample_distribution": "distributions/sample_distribution_summary.csv",
    "transformed_matrix": "matrices/transformed_matrix.csv",
    "detection_matrix": "matrices/detection_matrix.csv",
    "pairwise_correlations": "replicates/pairwise_correlations.csv",
    "pearson_matrix": "replicates/pearson_correlation_matrix.csv",
    "spearman_matrix": "replicates/spearman_correlation_matrix.csv",
    "shared_feature_matrix": "replicates/shared_feature_count_matrix.csv",
    "within_condition_correlations": "replicates/within_condition_correlation_summary.csv",
    "detection_overlap": "replicates/detection_overlap.csv",
    "jaccard_matrix": "replicates/jaccard_matrix.csv",
    "sample_distance_matrix": "replicates/sample_distance_matrix.csv",
    "feature_condition_cv": "variability/feature_condition_cv.csv",
    "condition_cv_summary": "variability/condition_cv_summary.csv",
    "pca_scores": "pca/pca_scores.csv",
    "pca_loadings": "pca/pca_loadings.csv",
    "pca_variance": "pca/pca_variance.csv",
    "pca_summary": "pca/pca_summary.csv",
    "sample_diagnostics": "diagnostics/sample_diagnostics.csv",
    "warnings": "diagnostics/warnings.csv",
    "summary": "summary.csv",
}
PLOTS = ("detected_features", "missing_fraction", "sample_distributions", "sample_boxplots",
    "pearson_heatmap", "spearman_heatmap", "shared_feature_heatmap", "jaccard_heatmap",
    "pca", "condition_cv", "feature_detection_frequency")


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""): digest.update(chunk)
    return digest.hexdigest()


def _configured_samples(project, selected_columns):
    if not selected_columns: raise NoQuantitativeColumnsError("Select quantitative columns for this QC run.")
    if len(selected_columns) < 2: raise TooFewSamplesError("At least two quantitative samples are required.",
        {"selected_sample_count": len(selected_columns)})
    if len(set(selected_columns)) != len(selected_columns):
        raise ProteomicsQCError("The same quantitative column was selected more than once.")
    configured = {item["column"]: item for item in experimental_design(project)["quantification_columns"]}
    missing = [name for name in selected_columns if name not in configured]
    if missing: raise NoQuantitativeColumnsError("Selected columns are not configured as Quantification.", {"columns": missing})
    types = {str(configured[name]["quantification_type"]) for name in selected_columns}
    if len(types) != 1 or "unknown" in types:
        raise MixedQuantificationTypesError("Select one explicitly configured quantification family.", {"types": sorted(types)})
    family = next(iter(types))
    if family not in DEFAULT_TRANSFORMATION:
        raise MixedQuantificationTypesError("Unsupported quantification family.", {"types": [family]})
    rows = []
    for index, name in enumerate(selected_columns, 1):
        item = configured[name]
        if not str(item["condition"]).strip():
            raise MissingConditionError(f"Condition metadata is missing for {name}.", {"sample": name})
        if not str(item["replicate"]).strip():
            raise MissingReplicateError(f"Replicate metadata is missing for {name}.", {"sample": name})
        rows.append({"sample_id": f"S{index:03d}", "column_name": name,
            "condition": str(item["condition"]), "replicate": str(item["replicate"]),
            "quantification_type": family})
    return pd.DataFrame(rows), family


def _validate_matrix(project, sample_metadata, family, transformation):
    processed = project.config.get("input", {}).get("processed_file")
    if not processed: raise ProteomicsQCError("Import a quantitative table before running QC.")
    path = Path(project.root) / processed
    if not path.is_file(): raise ProteomicsQCError("The processed input table is unavailable.")
    source = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = sample_metadata.column_name.tolist()
    missing = [name for name in columns if name not in source]
    if missing: raise NoQuantitativeColumnsError("Configured quantitative columns are absent from the input.", {"columns": missing})
    if source.empty: raise NoUsableQuantitativeDataError("The input table contains no features.")
    matrix = pd.DataFrame({"feature_id": [f"ROW{i:06d}" for i in range(1, len(source)+1)]})
    invalid = []
    negative = []
    usable = 0
    for sample in sample_metadata.itertuples(index=False):
        values = []
        for row_index, raw in enumerate(source[sample.column_name].astype(str), 1):
            stripped = raw.strip()
            if not stripped or stripped.casefold() in {"na", "nan", "null"}:
                values.append(""); continue
            try: number = float(stripped)
            except ValueError:
                invalid.append({"sample": sample.column_name, "source_row": row_index, "value": stripped})
                values.append(""); continue
            if math.isnan(number): values.append(""); continue
            if math.isinf(number): values.append("Inf" if number > 0 else "-Inf"); continue
            if number < 0 and (family in {"lfq_intensity", "raw_intensity", "spectral_count"} or transformation == "log2p1"):
                negative.append({"sample": sample.column_name, "source_row": row_index, "value": stripped})
            values.append(stripped)
            if transformation == "none" or number > 0 or (transformation == "log2p1" and number >= 0): usable += 1
        matrix[sample.sample_id] = values
    if invalid:
        raise InvalidNumericValueError("Text was found in declared quantitative columns.",
            {"invalid_value_count": len(invalid), "values": invalid})
    if negative:
        raise NegativeValueError("Negative values are incompatible with the selected quantification scale.",
            {"negative_value_count": len(negative), "values": negative})
    if not usable: raise NoUsableQuantitativeDataError("No usable quantitative values remain for this transformation.")
    primary = experimental_design(project)["primary_identifier_column"]
    identifiers = source[primary].astype(str) if primary in source else pd.Series([""] * len(source))
    feature = pd.DataFrame({"source_row": range(1, len(source)+1), "feature_id": matrix.feature_id,
        "display_identifier": identifiers,
        "original_identifier": identifiers})
    return matrix, feature


def prepare_qc_run(project, run_id, parameters: QCParameters):
    samples, family = _configured_samples(project, parameters.selected_columns)
    transformation = parameters.transformation or DEFAULT_TRANSFORMATION[family]
    if transformation not in TRANSFORMATIONS:
        raise UnsupportedTransformationError("Unsupported QC transformation.", {"transformation": transformation})
    if family in {"lfq_intensity", "raw_intensity"} and transformation == "log2p1":
        raise UnsupportedTransformationError("log2p1 is not permitted for linear intensity data.")
    zero_missing = parameters.zero_is_missing
    if zero_missing is None:
        if family not in DEFAULT_ZERO_MISSING:
            raise ProteomicsQCError("Specify zero_is_missing for other quantitative measurements.")
        zero_missing = DEFAULT_ZERO_MISSING[family]
    if not isinstance(zero_missing, bool): raise ProteomicsQCError("zero_is_missing must be Boolean.")
    matrix, features = _validate_matrix(project, samples, family, transformation)
    run = Path(project.root) / "analyses/Proteomics_QC/runs" / run_id
    provenance = Path(project.root) / "scripts/runs" / f"{run_id}_proteomics_qc"
    if run.exists() or provenance.exists(): raise ProteomicsQCError("This QC run ID already exists.")
    (run / "input").mkdir(parents=True)
    provenance.mkdir(parents=True)
    matrix.to_csv(run / "input/quantitative_matrix.csv", index=False)
    samples.to_csv(run / "input/sample_metadata.csv", index=False)
    features.to_csv(run / "input/feature_metadata.csv", index=False)
    hashes = {name: _sha(run / "input" / name) for name in
        ("quantitative_matrix.csv", "sample_metadata.csv", "feature_metadata.csv")}
    metadata = {"run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(),
        "quantification_type": family, "transformation": transformation,
        "zero_is_missing": zero_missing, "sample_count": len(samples), "feature_count": len(features),
        "selected_columns": list(parameters.selected_columns), "conditions": sorted(samples.condition.unique().tolist()),
        "sample_mapping": samples.to_dict("records"), "input_sha256": hashes,
        "feature_unit": "one original experimental row", "pca_center": True, "pca_scale": False,
        "pca_missingness_policy": "complete-case features only; no imputation",
        "network_access": False, "automatic_sample_exclusion": False}
    for path in (run / "input/parameters.json", provenance / "parameters.json"):
        path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    shutil.copy2(run / "input/sample_metadata.csv", provenance / "sample_metadata.csv")
    scripts = Path(__file__).resolve().parents[3] / "r_scripts"
    for source in (scripts / "11_proteomics_qc.R", scripts / "lib/proteomics_qc.R"):
        shutil.copy2(source, provenance / source.name)
    return run, provenance


def list_runs(project):
    folder = Path(project.root) / "analyses/Proteomics_QC/runs"
    return sorted((path.name for path in folder.iterdir() if path.is_dir()), reverse=True) if folder.is_dir() else []


def load_run(project, run_id):
    run = Path(project.root) / "analyses/Proteomics_QC/runs" / run_id
    required = {name: run / relative for name, relative in TABLES.items()}
    for name in ("quantitative_matrix.csv", "sample_metadata.csv", "feature_metadata.csv", "parameters.json"):
        required[f"input_{name}"] = run / "input" / name
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing: raise MissingQCOutputError("Missing required QC output.", {"paths": missing})
    try:
        metadata = json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))
        tables = {name: pd.read_csv(path, dtype=str, keep_default_na=False)
                  for name, path in required.items() if not name.startswith("input_")}
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise MissingQCOutputError(f"Invalid QC run output: {error}") from error
    return {"run_root": run, "metadata": metadata, "tables": tables,
        "workbook": (run / "Proteomics_QC.xlsx") if (run / "Proteomics_QC.xlsx").is_file() else None,
        "plots": tuple(sorted((run / "plots").glob("*.png")))}


def run_proteomics_qc(project, runtime: RRuntime, *, run_id, parameters: QCParameters,
                      script=None, timeout=300):
    run, provenance = prepare_qc_run(project, run_id, parameters)
    entry = script or Path(__file__).resolve().parents[3] / "r_scripts/11_proteomics_qc.R"
    try: result = runtime.run(entry, "--run", str(run), "--provenance", str(provenance), timeout=timeout)
    except RuntimeError as error: raise QCRExecutionError(str(error)) from error
    if result.returncode: raise QCRExecutionError((result.stderr or result.stdout or "Rscript failed.").strip())
    outputs = load_run(project, run_id)
    latest = Path(project.root) / "analyses/Proteomics_QC"
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "Proteomics_QC.xlsx", latest / "Proteomics_QC.xlsx")
    return outputs
