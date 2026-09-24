"""Offline limma analysis of an immutable Differential Preparation run."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import differential_preparation as preparation
from .r_runtime import RRuntime


class DifferentialAnalysisError(RuntimeError):
    code = "differential_analysis_error"
    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


class NoReadyPreparationError(DifferentialAnalysisError): code = "no_ready_preparation_run"
class PreparationRunIncompleteError(DifferentialAnalysisError): code = "preparation_run_incomplete"
class PreparedMatrixHashMismatchError(DifferentialAnalysisError): code = "prepared_matrix_hash_mismatch"
class SampleMetadataMismatchError(DifferentialAnalysisError): code = "sample_metadata_mismatch"
class ConditionMetadataMismatchError(DifferentialAnalysisError): code = "condition_metadata_mismatch"
class PreparedScaleRequiredError(DifferentialAnalysisError): code = "prepared_scale_required"
class UnsupportedPreparedScaleError(DifferentialAnalysisError): code = "unsupported_prepared_scale"
class DesignNotFullRankError(DifferentialAnalysisError): code = "design_not_full_rank"
class NoPreparedFeaturesError(DifferentialAnalysisError): code = "no_prepared_features"
class NoEstimableTestsError(DifferentialAnalysisError): code = "no_estimable_differential_tests"
class LimmaUnavailableError(DifferentialAnalysisError): code = "limma_package_unavailable"
class StatmodUnavailableError(DifferentialAnalysisError): code = "statmod_package_unavailable"
class DifferentialRExecutionError(DifferentialAnalysisError): code = "r_execution_failure"
class MissingDifferentialOutputError(DifferentialAnalysisError): code = "missing_required_output"


@dataclass(frozen=True)
class DifferentialParameters:
    preparation_run_id: str
    prepared_scale: str | None = None
    ebayes_trend: bool = False
    ebayes_robust: bool = False
    fdr_threshold: float = 0.05
    minimum_absolute_effect: float = 0.0
    top_n: int = 20
    volcano_label_top_n: int = 20


TABLES = {
    "summary": "summary.csv",
    "design": "design/design_matrix.csv",
    "all_results": "results/differential_results_all.csv",
    "tested_results": "results/differential_results_tested.csv",
    "untested_results": "results/differential_results_untested.csv",
    "fdr_significant_results": "results/differential_results_fdr_significant.csv",
    "significant_results": "results/differential_results_significant.csv",
    "ranked_results": "results/differential_results_ranked.csv",
    "qualitative_candidates": "qualitative/qualitative_detection_candidates.csv",
    "preparation_metadata": "audit/preparation_run_metadata.csv",
    "model_audit": "audit/model_feature_audit.csv",
    "ebayes_audit": "audit/ebayes_summary.csv",
}
PLOTS = ("volcano", "ma", "mean_a_vs_b", "p_value_histogram", "fdr_histogram",
    "effect_distribution", "residual_df_distribution", "top_effects")


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parent(project, run_id):
    if not preparation.list_runs(project):
        raise NoReadyPreparationError("No Ready Differential Preparation run is available.")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", str(run_id)):
        raise PreparationRunIncompleteError("Invalid preparation run ID.")
    root = Path(project.root) / "analyses/Differential/Preparation/runs" / run_id
    parent_manifest = root / "input/output_hashes.json"
    prepared = root / "matrices/prepared_matrix.csv"
    if parent_manifest.is_file() and prepared.is_file():
        try:
            expected = json.loads(parent_manifest.read_text(encoding="utf-8")).get("prepared_matrix.csv")
        except (OSError, ValueError) as error:
            raise PreparationRunIncompleteError("Preparation output integrity manifest is invalid.") from error
        if expected != _sha(prepared):
            raise PreparedMatrixHashMismatchError("Prepared matrix hash mismatch.")
    try:
        output = preparation.load_run(project, run_id)
    except preparation.MissingPreparationOutputError as error:
        raise PreparationRunIncompleteError(f"Differential Preparation run is incomplete: {error}") from error
    run = output["run_root"]
    manifest = run / "input/output_hashes.json"
    if not manifest.is_file():
        raise PreparationRunIncompleteError("Preparation run has no output integrity manifest.")
    try:
        hashes = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise PreparationRunIncompleteError("Preparation output integrity manifest is invalid.") from error
    for name, relative in (("prepared_matrix.csv", "matrices/prepared_matrix.csv"),
        ("contrast_metadata.csv", "design/contrast_metadata.csv"),
        ("feature_eligibility.csv", "eligibility/feature_eligibility.csv"),
        ("qualitative_detection_candidates.csv", "eligibility/qualitative_detection_candidates.csv")):
        path = run / relative
        if hashes.get(name) != _sha(path):
            raise PreparedMatrixHashMismatchError("Preparation output hash mismatch.", {"artifact": name})
    inputs = output["metadata"].get("input_sha256", {})
    for name in ("sample_metadata.csv", "feature_metadata.csv", "original_quantitative_matrix.csv"):
        path = run / "input" / name
        if inputs.get(name) != _sha(path):
            raise PreparedMatrixHashMismatchError("Preparation input hash mismatch.", {"artifact": name})
    return output, hashes


def list_preparation_runs(project):
    ready = []
    for run_id in preparation.list_runs(project):
        try:
            _parent(project, run_id)
            ready.append(run_id)
        except DifferentialAnalysisError:
            pass
    return ready


def _validate_parent(output):
    meta = output["metadata"]
    tables = output["tables"]
    sample = pd.read_csv(output["run_root"] / "input/sample_metadata.csv", dtype=str, keep_default_na=False)
    feature = pd.read_csv(output["run_root"] / "input/feature_metadata.csv", dtype=str, keep_default_na=False)
    matrix = tables["prepared_matrix"]
    if matrix.empty:
        raise NoPreparedFeaturesError("The parent prepared matrix has no features.")
    if sample.sample_id.duplicated().any() or feature.feature_id.duplicated().any() or matrix.feature_id.duplicated().any():
        raise SampleMetadataMismatchError("Duplicated sample or feature IDs in the parent run.")
    expected = sample.sample_id.tolist()
    actual = [column for column in matrix if column.startswith("S") and column[1:].isdigit()]
    if actual != expected:
        raise SampleMetadataMismatchError("Prepared matrix samples do not match frozen sample metadata.")
    conditions = sample.condition.unique().tolist()
    if conditions != [meta["condition_A"], meta["condition_B"]] or len(conditions) != 2:
        raise ConditionMetadataMismatchError("Frozen preparation condition order is inconsistent.")
    contrast = tables["contrast_metadata"].iloc[0]
    if contrast.condition_A != meta["condition_A"] or contrast.condition_B != meta["condition_B"]:
        raise ConditionMetadataMismatchError("Frozen contrast metadata differs from preparation parameters.")
    eligibility = tables["feature_eligibility"]
    eligible_ids = eligibility.loc[eligibility.continuous_eligible == "TRUE", "feature_id"].tolist()
    if matrix.feature_id.tolist() != eligible_ids:
        raise PreparationRunIncompleteError("Prepared matrix feature IDs do not match frozen eligibility.")
    if not set(matrix.feature_id).issubset(set(feature.feature_id)):
        raise PreparationRunIncompleteError("Prepared matrix feature metadata is incomplete.")
    return sample, feature


def prepare_run(project, run_id, parameters: DifferentialParameters):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", str(run_id)):
        raise DifferentialAnalysisError("Run ID must be a safe filename component.")
    if not isinstance(parameters.ebayes_trend, bool) or not isinstance(parameters.ebayes_robust, bool):
        raise DifferentialAnalysisError("eBayes trend and robust options must be Boolean.")
    if not 0 < parameters.fdr_threshold <= 1 or parameters.minimum_absolute_effect < 0:
        raise DifferentialAnalysisError("Invalid FDR or minimum effect threshold.")
    if parameters.top_n < 1 or parameters.volcano_label_top_n < 0:
        raise DifferentialAnalysisError("Invalid plotting limits.")
    parent, parent_hashes = _parent(project, parameters.preparation_run_id)
    samples, _ = _validate_parent(parent)
    transformation = parent["metadata"]["transformation"]
    if transformation == "log2_positive":
        if parameters.prepared_scale not in (None, "log2"):
            raise UnsupportedPreparedScaleError("log2-positive preparation must be analysed as log2 scale.")
        scale = "log2"
    else:
        if parameters.prepared_scale is None:
            raise PreparedScaleRequiredError("Declare prepared_scale as log2 or continuous for untransformed data.")
        if parameters.prepared_scale not in ("log2", "continuous"):
            raise UnsupportedPreparedScaleError("Unsupported prepared scale.")
        scale = parameters.prepared_scale
    run = Path(project.root) / "analyses/Differential/Statistics/runs" / run_id
    provenance = Path(project.root) / "scripts/runs" / f"{run_id}_differential_statistics"
    if run.exists() or provenance.exists():
        raise DifferentialAnalysisError("This differential run ID already exists.")
    (run / "input").mkdir(parents=True)
    provenance.mkdir(parents=True)
    copies = {
        "prepared_matrix.csv": "matrices/prepared_matrix.csv",
        "sample_metadata.csv": "input/sample_metadata.csv",
        "feature_metadata.csv": "input/feature_metadata.csv",
        "contrast_metadata.csv": "design/contrast_metadata.csv",
        "preparation_summary.csv": "summary.csv",
        "preparation_parameters.json": "input/parameters.json",
        "feature_eligibility.csv": "eligibility/feature_eligibility.csv",
        "qualitative_detection_candidates.csv": "eligibility/qualitative_detection_candidates.csv",
        "feature_imputation_summary.csv": "imputation/feature_imputation_summary.csv",
        "imputation_mask.csv": "imputation/imputation_mask.csv",
    }
    for destination, source in copies.items():
        shutil.copy2(parent["run_root"] / source, run / "input" / destination)
    parent_meta = parent["metadata"]
    metadata = {"run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(),
        "parent_preparation_run_id": parameters.preparation_run_id,
        "preparation_run_id": parameters.preparation_run_id,
        "preparation_run_path": str(parent["run_root"]),
        "preparation_parameters_hash": _sha(run / "input/preparation_parameters.json"),
        "prepared_matrix_hash": _sha(run / "input/prepared_matrix.csv"),
        "contrast_metadata_hash": _sha(run / "input/contrast_metadata.csv"),
        "parent_output_hashes": parent_hashes,
        "condition_A": parent_meta["condition_A"], "condition_B": parent_meta["condition_B"],
        "comparison_direction": "Condition A - Condition B",
        "prepared_scale": scale, "effect_scale": scale,
        "preparation_transformation": parent_meta["transformation"],
        "preparation_normalization": parent_meta["normalization"],
        "preparation_imputation": parent_meta["imputation_method"],
        "preparation_imputation_seed": parent_meta["imputation_seed"],
        "preparation_zero_is_missing": parent_meta["zero_is_missing"],
        "preparation_minimum_observed": parent_meta["minimum_observed_per_condition"],
        "sample_count_A": int((samples.condition == parent_meta["condition_A"]).sum()),
        "sample_count_B": int((samples.condition == parent_meta["condition_B"]).sum()),
        "ebayes_trend": parameters.ebayes_trend, "ebayes_robust": parameters.ebayes_robust,
        "fdr_method": "BH", "fdr_threshold": parameters.fdr_threshold,
        "minimum_absolute_effect": parameters.minimum_absolute_effect,
        "confidence_level": 0.95, "top_n": parameters.top_n,
        "volcano_label_top_n": parameters.volcano_label_top_n,
        "imputation_caveat": "Imputed values are treated as quantitative inputs to the fitted linear model. The differential model does not propagate imputation uncertainty.",
        "feature_unit": "one original experimental row", "network_access": False}
    for destination in (run / "input/parameters.json", provenance / "parameters.json"):
        destination.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    shutil.copy2(run / "input/preparation_parameters.json", provenance / "preparation_parameters.json")
    scripts = Path(__file__).resolve().parents[3] / "r_scripts"
    for source in (scripts / "13_differential_analysis.R", scripts / "lib/differential_analysis.R"):
        shutil.copy2(source, provenance / source.name)
    return run, provenance


def list_runs(project):
    folder = Path(project.root) / "analyses/Differential/Statistics/runs"
    return sorted((path.name for path in folder.iterdir() if path.is_dir()), reverse=True) if folder.is_dir() else []


def load_run(project, run_id):
    run = Path(project.root) / "analyses/Differential/Statistics/runs" / run_id
    required = {name: run / relative for name, relative in TABLES.items()}
    for name in ("prepared_matrix.csv", "sample_metadata.csv", "feature_metadata.csv",
        "contrast_metadata.csv", "preparation_summary.csv", "preparation_parameters.json",
        "feature_eligibility.csv", "qualitative_detection_candidates.csv", "parameters.json"):
        required[f"input_{name}"] = run / "input" / name
    required["workbook"] = run / "Differential_analysis.xlsx"
    for plot in PLOTS:
        for suffix in (".png", ".pdf"):
            required[f"plot_{plot}{suffix}"] = run / "plots" / f"{plot}{suffix}"
    parameters_path = run / "input/parameters.json"
    if parameters_path.is_file():
        try:
            if json.loads(parameters_path.read_text(encoding="utf-8")).get("preparation_imputation") != "none":
                for suffix in (".png", ".pdf"):
                    required[f"plot_imputation_diagnostic{suffix}"] = run / "plots" / f"imputation_diagnostic{suffix}"
        except (OSError, ValueError) as error:
            raise MissingDifferentialOutputError("Invalid run parameters.") from error
    absent = [str(path) for path in required.values() if not path.is_file()]
    if absent:
        raise MissingDifferentialOutputError("Differential run is missing required artifacts.", {"paths": absent})
    try:
        metadata = json.loads((run / "input/parameters.json").read_text(encoding="utf-8"))
        tables = {name: pd.read_csv(path, dtype=str, keep_default_na=False)
            for name, path in required.items() if name in TABLES}
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise MissingDifferentialOutputError(f"Could not load differential run: {error}") from error
    if metadata["prepared_matrix_hash"] != _sha(run / "input/prepared_matrix.csv"):
        raise PreparedMatrixHashMismatchError("Frozen differential input hash mismatch.")
    return {"run_root": run, "metadata": metadata, "tables": tables,
        "plots": tuple(sorted((run / "plots").glob("*.png"))),
        "workbook": run / "Differential_analysis.xlsx"}


def run_differential_analysis(project, runtime: RRuntime, *, run_id,
                              parameters: DifferentialParameters, script=None, timeout=300):
    run, provenance = prepare_run(project, run_id, parameters)
    entry = script or Path(__file__).resolve().parents[3] / "r_scripts/13_differential_analysis.R"
    try:
        result = runtime.run(entry, "--run", str(run), "--provenance", str(provenance), timeout=timeout)
    except RuntimeError as error:
        raise DifferentialRExecutionError(str(error)) from error
    if result.returncode:
        message = (result.stderr or result.stdout or "Rscript failed.").strip()
        if "LIMMA_UNAVAILABLE" in message:
            raise LimmaUnavailableError("The limma R package is unavailable.")
        if "STATMOD_UNAVAILABLE" in message:
            raise StatmodUnavailableError("The statmod R package is required for robust empirical Bayes.")
        if "DESIGN_NOT_FULL_RANK" in message:
            raise DesignNotFullRankError("The two-condition design is not full rank.")
        if "NO_ESTIMABLE_TESTS" in message:
            raise NoEstimableTestsError("No features produced an estimable differential test.")
        raise DifferentialRExecutionError(message)
    outputs = load_run(project, run_id)
    latest = Path(project.root) / "analyses/Differential/Statistics"
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "Differential_analysis.xlsx", latest / "Differential_analysis.xlsx")
    return outputs
