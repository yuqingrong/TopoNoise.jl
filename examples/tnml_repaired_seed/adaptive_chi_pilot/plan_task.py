#!/usr/bin/env python3
"""Read-only phase-1 task planner; numerical execution is deliberately disabled."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

try:  # Package import.
    from . import manifest
except ImportError:  # Direct CLI invocation.
    import manifest


def plan_task(root: str | Path, stage: str, task_id: int) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    document = manifest.load_manifest(root_path)
    task = manifest.task_for_id(document, stage, task_id)
    output = (root_path / task["output"]).resolve()
    if output == root_path or root_path not in output.parents:
        raise ValueError("task output resolves outside the study root")
    return {
        "status": "dry_run",
        "execution": "disabled_phase1",
        "manifest_id": document["manifest_id"],
        "stage": stage,
        "task_id": task_id,
        "task": task,
        "output_directory": str(output),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Print one immutable adaptive-chi task without executing it."
    )
    parser.add_argument("--root", required=True)
    parser.add_argument("--stage", required=True, choices=("l3_exact", "l5_feasibility"))
    parser.add_argument("--task-id", required=True, type=int)
    arguments = parser.parse_args(argv)
    payload = plan_task(arguments.root, arguments.stage, arguments.task_id)
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
