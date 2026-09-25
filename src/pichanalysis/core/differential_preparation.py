"""Freeze and load offline quantitative matrices for later differential modelling.

No model fitting or differential statistics are performed here.
"""
from __future__ import annotations
from .resources import r_script, r_scripts_dir

import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .column_mapping import experimental_design
from .r_runtime import RRuntime


class DifferentialPreparationError(RuntimeError):
    code = "differential_preparation_error"
    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


class TooFewConditionsError(DifferentialPreparationError): code = "fewer_than_two_conditions"
class SameConditionError(DifferentialPreparationError): code = "same_condition_selected_twice"
class TooFewSamplesError(DifferentialPreparationError): code = "fewer_than_two_samples_in_condition"
class MixedQuantificationTypesError(DifferentialPreparationError): code = "mixed_quantification_types"
class SpectralCountsUnsupportedError(DifferentialPreparationError): code = "spectral_counts_not_supported"
class OtherContinuousUnconfirmedError(DifferentialPreparationError): code = "other_continuous_not_confirmed"
class MissingSampleMetadataError(DifferentialPreparationError): code = "missing_sample_metadata"
class InvalidNumericValueError(DifferentialPreparationError): code = "invalid_numeric_value"
class NegativeValueError(DifferentialPreparationError): code = "negative_value_incompatible"
class ZeroLog2IncompatibleError(DifferentialPreparationError): code = "zero_log2_incompatible"
class UnsupportedTransformationError(DifferentialPreparationError): code = "unsupported_transformation"
class UnsupportedNormalizationError(DifferentialPreparationError): code = "unsupported_normalization"
class UnsupportedImputationError(DifferentialPreparationError): code = "unsupported_imputation_method"
class ImputationPackageUnavailableError(DifferentialPreparationError): code = "imputation_package_unavailable"
class NoEligibleFeaturesError(DifferentialPreparationError): code = "no_continuous_eligible_features"
class PreparationRExecutionError(DifferentialPreparationError): code = "r_execution_failure"
class MissingPreparationOutputError(DifferentialPreparationError): code = "missing_required_output"


@dataclass(frozen=True)
class PreparationParameters:
    condition_a: str
    condition_b: str
    selected_columns: tuple[str, ...] = ()
    other_continuous_confirmed: bool = False
    transformation: str | None = None
    zero_is_missing: bool | None = None
    normalization: str = "none"
    minimum_observed_per_condition: int = 2
    imputation_method: str = "none"
    imputation_seed: int = 12345
    minprob_q: float = 0.01
    imputation_sigma: float = 1.0
    imputation_margin: int = 2
    knn_k: int = 10


TABLES = {
    "summary": "summary.csv",
    "original_linear_matrix": "matrices/original_linear_matrix.csv",
    "transformed_matrix": "matrices/transformed_matrix.csv",
    "normalized_matrix": "matrices/normalized_matrix.csv",
    "prepared_matrix": "matrices/prepared_matrix.csv",
    "feature_missingness_patterns": "eligibility/feature_missingness_patterns.csv",
    "feature_eligibility": "eligibility/feature_eligibility.csv",
    "qualitative_detection_candidates": "eligibility/qualitative_detection_candidates.csv",
    "excluded_features": "eligibility/excluded_features.csv",
    "sample_normalization": "normalization/sample_normalization.csv",
    "imputation_mask": "imputation/imputation_mask.csv",
    "imputation_method_parameters": "imputation/method_parameters.csv",
    "imputed_cells": "imputation/imputed_cells.csv",
    "feature_imputation_summary": "imputation/feature_imputation_summary.csv",
    "unresolved_missing_values": "imputation/unresolved_missing_values.csv",
    "contrast_metadata": "design/contrast_metadata.csv",
    "design_preview": "design/design_preview.csv",
}
PLOTS = ("missingness_patterns", "observed_replicates", "eligibility_counts",
    "before_normalization", "after_normalization", "before_imputation_missingness",
    "imputed_per_sample", "before_after_imputation", "qualitative_candidates")


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _select_samples(project, parameters):
    entries = experimental_design(project)["quantification_columns"]
    conditions = {str(entry["condition"]).strip() for entry in entries if entry["condition"]}
    if len(conditions) < 2:
        raise TooFewConditionsError("At least two configured quantitative conditions are required.")
    if parameters.condition_a == parameters.condition_b:
        raise SameConditionError("Condition A and Condition B must differ.")
    if parameters.condition_a not in conditions or parameters.condition_b not in conditions:
        raise TooFewConditionsError("Both selected conditions must be available in the project.")
    relevant = [entry for entry in entries if entry["condition"] in (parameters.condition_a, parameters.condition_b)]
    by_name = {entry["column"]: entry for entry in relevant}
    selected = parameters.selected_columns or tuple(entry["column"] for entry in relevant)
    if len(selected) != len(set(selected)) or any(name not in by_name for name in selected):
        raise MissingSampleMetadataError("Selected samples must be distinct, configured quantitative columns in the selected conditions.")
    chosen = [by_name[name] for name in selected]
    for entry in chosen:
        if not str(entry["condition"]).strip() or not str(entry["replicate"]).strip():
            raise MissingSampleMetadataError("Selected samples require condition and replicate metadata.", {"sample": entry["column"]})
    counts = {condition: sum(entry["condition"] == condition for entry in chosen)
        for condition in (parameters.condition_a, parameters.condition_b)}
    if min(counts.values()) < 2:
        raise TooFewSamplesError("Select at least two samples in each selected condition.", {"sample_counts": counts})
    families = {str(entry["quantification_type"]) for entry in chosen}
    if len(families) != 1 or "unknown" in families:
        raise MixedQuantificationTypesError("Selected samples use incompatible quantification types.", {"types": sorted(families)})
    family = next(iter(families))
    if family == "spectral_count":
        raise SpectralCountsUnsupportedError("Spectral-count differential analysis requires a count-based model and is not supported by the current intensity-based differential pipeline.")
    if family == "other_quantitative" and not parameters.other_continuous_confirmed:
        raise OtherContinuousUnconfirmedError("Confirm that Other quantitative values are continuous before using this pipeline.")
    if family not in {"lfq_intensity", "raw_intensity", "other_quantitative"}:
        raise MixedQuantificationTypesError("Unsupported quantitative family.", {"types": sorted(families)})
    ordered = sorted(chosen, key=lambda entry: 0 if entry["condition"] == parameters.condition_a else 1)
    samples = pd.DataFrame([{"sample_id": f"S{index:03d}", "column_name": entry["column"],
        "condition": entry["condition"], "replicate": entry["replicate"],
        "quantification_type": family} for index, entry in enumerate(ordered, 1)])
    return samples, family


def _options(parameters, family):
    transformation = parameters.transformation or ("none" if family == "other_quantitative" else "log2_positive")
    if transformation not in {"none", "log2_positive"}:
        raise UnsupportedTransformationError("Choose log2_positive or none.")
    if parameters.normalization not in {"none", "median_center"}:
        raise UnsupportedNormalizationError("Choose none or median_center normalization.")
    if parameters.imputation_method not in {"none", "MinProb", "QRILC", "KNN"}:
        raise UnsupportedImputationError("Choose none, MinProb, QRILC, or KNN imputation.")
    if not isinstance(parameters.minimum_observed_per_condition, int) or parameters.minimum_observed_per_condition < 1:
        raise DifferentialPreparationError("minimum_observed_per_condition must be a positive integer.")
    if not isinstance(parameters.imputation_seed, int) or parameters.imputation_seed < 0:
        raise DifferentialPreparationError("imputation_seed must be a nonnegative integer.")
    if not 0 < parameters.minprob_q < 1 or parameters.imputation_sigma <= 0 or parameters.imputation_margin not in (1, 2) or parameters.knn_k < 1:
        raise DifferentialPreparationError("Invalid imputation parameters.")
    zero_missing = parameters.zero_is_missing
    if zero_missing is None:
        if family == "other_quantitative":
            raise DifferentialPreparationError("Set zero_is_missing explicitly for Other quantitative values.")
        zero_missing = True
    if not isinstance(zero_missing, bool):
        raise DifferentialPreparationError("zero_is_missing must be Boolean.")
    return transformation, zero_missing


def _freeze_matrix(project, samples, transformation, zero_missing):
    relative = project.config.get("input", {}).get("processed_file")
    if not relative:
        raise DifferentialPreparationError("Import a quantitative table first.")
    source = pd.read_csv(Path(project.root) / relative, dtype=str, keep_default_na=False)
    if source.empty:
        raise NoEligibleFeaturesError("The input table has no features.")
    missing = [name for name in samples.column_name if name not in source]
    if missing:
        raise MissingSampleMetadataError("Configured quantitative columns are absent from the input.", {"columns": missing})
    invalid, negative, zeros = [], [], []
    value_counts = {"blank": 0, "NA": 0, "NaN": 0, "null": 0, "Inf": 0, "-Inf": 0, "zero": 0}
    matrix = pd.DataFrame({"feature_id": [f"ROW{i:06d}" for i in range(1, len(source) + 1)]})
    for sample in samples.itertuples(index=False):
        values = []
        for row, raw in enumerate(source[sample.column_name].astype(str), 1):
            value = raw.strip()
            if not value or value.casefold() in {"na", "nan", "null"}:
                key = "blank" if not value else {"na": "NA", "nan": "NaN", "null": "null"}[value.casefold()]
                value_counts[key] += 1
                values.append("")
                continue
            try:
                number = float(value)
            except ValueError:
                invalid.append({"sample": sample.column_name, "source_row": row, "value": value})
                values.append(value)
                continue
            if number == 0:
                value_counts["zero"] += 1
            if not math.isfinite(number):
                value_counts["Inf" if number > 0 else "-Inf"] += 1
            elif number < 0 and transformation == "log2_positive":
                negative.append({"sample": sample.column_name, "source_row": row, "value": value})
            elif number == 0 and transformation == "log2_positive" and not zero_missing:
                zeros.append({"sample": sample.column_name, "source_row": row})
            values.append(value)
        matrix[sample.sample_id] = values
    if invalid:
        raise InvalidNumericValueError("Text was found in a quantitative column.",
            {"invalid_value_count": len(invalid), "values": invalid})
    if negative:
        raise NegativeValueError("Negative values are incompatible with log2_positive.", {"values": negative})
    if zeros:
        raise ZeroLog2IncompatibleError("Zero is a legitimate quantitative value under this policy but cannot be log2-transformed without a pseudocount.", {"values": zeros})
    primary = experimental_design(project)["primary_identifier_column"]
    identifiers = source[primary].astype(str) if primary in source else pd.Series([""] * len(source))
    features = pd.DataFrame({"feature_id": matrix.feature_id, "source_row": range(1, len(source) + 1),
        "display_identifier": identifiers, "original_identifier": identifiers})
    for identifier_type, output_name in (("gene_symbol", "gene_symbol"),
                                         ("uniprot", "uniprot_accession")):
        configured = next((name for name, config in project.config.get("columns", {}).items()
            if name in source and config.get("role") == "identifier"
            and config.get("identifier_type") == identifier_type), None)
        features[output_name] = source[configured].astype(str) if configured else ""
    return matrix, features, value_counts


def prepare_run(project, run_id, parameters: PreparationParameters):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", str(run_id)):
        raise DifferentialPreparationError("Run ID must be a safe filename component.")
    samples, family = _select_samples(project, parameters)
    transformation, zero_missing = _options(parameters, family)
    matrix, features, value_counts = _freeze_matrix(project, samples, transformation, zero_missing)
    run = Path(project.root) / "analyses/Differential/Preparation/runs" / run_id
    provenance = Path(project.root) / "scripts/runs" / f"{run_id}_differential_preparation"
    if run.exists() or provenance.exists():
        raise DifferentialPreparationError("This preparation run ID already exists.")
    (run / "input").mkdir(parents=True)
    provenance.mkdir(parents=True)
    paths = {"original_quantitative_matrix.csv": matrix, "sample_metadata.csv": samples,
        "feature_metadata.csv": features}
    for name, frame in paths.items():
        frame.to_csv(run / "input" / name, index=False)
    hashes = {name: _sha(run / "input" / name) for name in paths}
    metadata = {"run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(),
        "condition_A": parameters.condition_a, "condition_B": parameters.condition_b,
        "comparison_direction": "Condition A - Condition B on log scale; Condition A / Condition B on linear scale",
        "sample_count_A": int((samples.condition == parameters.condition_a).sum()),
        "sample_count_B": int((samples.condition == parameters.condition_b).sum()),
        "selected_columns": samples.column_name.tolist(), "sample_mapping": samples.to_dict("records"),
        "quantification_type": family, "other_continuous_confirmed": parameters.other_continuous_confirmed,
        "transformation": transformation, "zero_is_missing": zero_missing,
        "normalization": parameters.normalization,
        "minimum_observed_per_condition": parameters.minimum_observed_per_condition,
        "imputation_method": parameters.imputation_method, "imputation_seed": parameters.imputation_seed,
        "imputation_scope": "partial missing values in continuous-analysis eligible features only",
        "imputation_parameters": {"q": parameters.minprob_q, "sigma": parameters.imputation_sigma,
            "MARGIN": parameters.imputation_margin, "knn_k": parameters.knn_k},
        "feature_unit": "one original experimental row", "input_sha256": hashes,
        "input_value_counts": value_counts, "network_access": False,
        "automatic_sample_exclusion": False, "missing_mechanism_inferred": False}
    for destination in (run / "input/parameters.json", provenance / "parameters.json"):
        destination.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    shutil.copy2(run / "input/sample_metadata.csv", provenance / "sample_metadata.csv")
    scripts = r_scripts_dir()
    for source in (scripts / "12_differential_preparation.R", scripts / "lib/differential_preparation.R"):
        shutil.copy2(source, provenance / source.name)
    return run, provenance


def list_runs(project):
    root = Path(project.root) / "analyses/Differential/Preparation/runs"
    return sorted((path.name for path in root.iterdir() if path.is_dir()), reverse=True) if root.is_dir() else []


def load_run(project, run_id):
    run = Path(project.root) / "analyses/Differential/Preparation/runs" / run_id
    required = {name: run / relative for name, relative in TABLES.items()}
    for name in ("original_quantitative_matrix.csv", "sample_metadata.csv", "feature_metadata.csv", "parameters.json"):
        required[f"input_{name}"] = run / "input" / name
    required["workbook"] = run / "Differential_preparation.xlsx"
    absent = [str(path) for path in required.values() if not path.is_file()]
    if absent:
        raise MissingPreparationOutputError("Required preparation artifacts are missing.", {"paths": absent})
    try:
        metadata = json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))
        tables = {name: pd.read_csv(path, dtype=str, keep_default_na=False)
            for name, path in required.items() if name in TABLES}
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise MissingPreparationOutputError(f"Could not load preparation run: {error}") from error
    return {"run_root": run, "metadata": metadata, "tables": tables,
        "plots": tuple(sorted((run / "plots").glob("*.png"))),
        "workbook": run / "Differential_preparation.xlsx"}


def run_differential_preparation(project, runtime: RRuntime, *, run_id,
                                 parameters: PreparationParameters, script=None, timeout=300):
    run, provenance = prepare_run(project, run_id, parameters)
    entry = script or r_script("12_differential_preparation.R")
    try:
        result = runtime.run(entry, "--run", str(run), "--provenance", str(provenance), timeout=timeout)
    except RuntimeError as error:
        raise PreparationRExecutionError(str(error)) from error
    if result.returncode:
        message = (result.stderr or result.stdout or "Rscript failed.").strip()
        if "NO_CONTINUOUS_ELIGIBLE" in message:
            raise NoEligibleFeaturesError("No features are eligible for continuous differential analysis. Qualitative detection candidates were identified separately.")
        if "IMPUTATION_PACKAGE_UNAVAILABLE" in message:
            raise ImputationPackageUnavailableError(message)
        raise PreparationRExecutionError(message)
    output_paths = {"prepared_matrix.csv": run / "matrices/prepared_matrix.csv",
        "contrast_metadata.csv": run / "design/contrast_metadata.csv",
        "feature_eligibility.csv": run / "eligibility/feature_eligibility.csv",
        "qualitative_detection_candidates.csv": run / "eligibility/qualitative_detection_candidates.csv"}
    output_hashes = {name: _sha(path) for name, path in output_paths.items()}
    for destination in (run / "input/output_hashes.json", provenance / "output_hashes.json"):
        destination.write_text(json.dumps(output_hashes, indent=2), encoding="utf-8")
    outputs = load_run(project, run_id)
    latest = Path(project.root) / "analyses/Differential/Preparation"
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "Differential_preparation.xlsx", latest / "Differential_preparation.xlsx")
    return outputs
