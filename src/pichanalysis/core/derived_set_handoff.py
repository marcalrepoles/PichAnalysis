"""Safe, offline handoff of frozen comparison entities to manual targets."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .experiment_comparison import IDENTITY_COLUMNS, ComparisonError, load_run

DERIVED_SETS = {
    "a_only": "A only", "b_only": "B only", "shared": "Shared",
    "significant_up_both": "Up in both", "significant_down_both": "Down in both",
    "opposite_direction": "Opposite direction",
    "significant_only_a": "Significant only in A",
    "significant_only_b": "Significant only in B",
}
DESTINATIONS = ("go", "kegg", "reactome", "mitocarta", "mtdna", "domains", "string", "complexes")


def prepare(project, comparison_run_id: str, derived_set: str, destination: str) -> dict:
    """Resolve a frozen set to current project source rows; never run analysis."""
    if derived_set not in DERIVED_SETS or destination not in DESTINATIONS:
        raise ComparisonError("Unsupported derived set or destination.")
    loaded = load_run(project, comparison_run_id)
    frame = loaded["tables"].get(derived_set)
    if frame is None or frame.empty:
        raise ComparisonError("The selected derived set is empty or unavailable.")
    identity = loaded["config"].get("comparison_type")
    column = IDENTITY_COLUMNS.get(identity)
    if column is None:
        raise ComparisonError("This comparison identity cannot be safely mapped to project source rows.")
    path = Path(project.root) / "mapping/tables/protein_catalog.csv"
    if not path.is_file():
        raise ComparisonError("The current project has no Mapping catalog for manual target selection.")
    catalog = pd.read_csv(path, dtype=str, keep_default_na=False)
    if column not in catalog or "source_row" not in catalog:
        raise ComparisonError("The current Mapping catalog lacks the required identifier or source row.")
    entities = set(frame.comparison_entity.astype(str).str.strip())
    if not entities or "" in entities:
        raise ComparisonError("The derived set has missing comparison identifiers.")
    matched = catalog[catalog[column].astype(str).str.strip().isin(entities)]
    found = set(matched[column].astype(str).str.strip())
    missing = sorted(entities - found)
    if missing:
        raise ComparisonError(f"Cannot hand off {len(missing)} entities absent from the current Mapping catalog.")
    try:
        rows = sorted({int(value) for value in matched.source_row})
    except ValueError as error:
        raise ComparisonError("The current Mapping catalog has invalid source rows.") from error
    if not rows:
        raise ComparisonError("No project source rows were resolved.")
    return {"source": "Experiment Comparison", "comparison_run_id": comparison_run_id,
        "derived_set": derived_set, "destination": destination, "comparison_identity": identity,
        "mapping_catalog": str(path), "resolved_entities": len(entities), "source_rows": rows,
        "prepared_at": datetime.now(timezone.utc).isoformat()}


def record(project, handoff: dict) -> Path:
    """Persist provenance only after the GUI has preloaded the manual target."""
    root = Path(project.root) / "analyses/experiment_comparison/runs" / handoff["comparison_run_id"] / "handoffs"
    root.mkdir(parents=True, exist_ok=True)
    import uuid
    identifier = uuid.uuid4().hex
    path = root / f"{identifier}.json"
    path.write_text(json.dumps(handoff, indent=2) + "\n", encoding="utf-8")
    destination = Path(project.root) / "analyses" / handoff["destination"] / "handoffs"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / f"{identifier}.json").write_text(
        json.dumps(handoff | {"status": "Awaiting user-confirmed analysis"}, indent=2) + "\n", encoding="utf-8")
    return path
