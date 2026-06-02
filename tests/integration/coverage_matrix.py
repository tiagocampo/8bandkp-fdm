#!/usr/bin/env python3
"""Validation Coverage Matrix generator.

Reads a YAML universe file declaring (observable, geometry, material) cells
and scans test source files for COVERAGE annotations. Cross-checks both
directions and produces a dual-layer markdown report.

Usage: coverage_matrix.py --source-dir <repo_root> [--output <file>]

Output: markdown report to stdout (or file), summary to stderr.
Exit: 1 if any required-tier cell has no test coverage, 0 otherwise.

Annotation format (in .sh and .py test files):
  # COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
  # COVERAGE: observable=m*_e geometry=bulk material=InAs
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import Any


UNIVERSE_FILE = "tests/integration/validation_universe.yml"
SCAN_PATTERNS = [
    "tests/integration/verify_*.py",
    "tests/integration/test_*.sh",
    "tests/integration/test_*.py",
]
ANNOTATION_RE = re.compile(r"^#\s*COVERAGE:\s*(.+)$", re.MULTILINE)
KV_RE = re.compile(r"(\w+)=([^\s]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validation Coverage Matrix")
    parser.add_argument("--source-dir", required=True, help="Repo root directory")
    parser.add_argument("--output", help="Output file path (default: stdout)")
    return parser.parse_args()


def load_universe(source_dir: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        import yaml  # noqa: F811
    except ImportError:
        print("Error: PyYAML is required. Install with: pip install pyyaml",
              file=sys.stderr)
        sys.exit(1)

    path = os.path.join(source_dir, UNIVERSE_FILE)
    with open(path) as f:
        data = yaml.safe_load(f)
    cells = []
    for i, c in enumerate(data.get("cells", [])):
        for required in ("observable", "geometry", "tier"):
            if required not in c:
                print(f"Error: cell {i} missing required field '{required}'", file=sys.stderr)
                sys.exit(1)
        cells.append({
            "observable": c["observable"],
            "geometry": c["geometry"],
            "material": c.get("material"),
            "tier": c["tier"],
            "reference": c.get("reference"),
        })
    return cells, data.get("metadata", {})


def scan_annotations(source_dir: str) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    for pattern in SCAN_PATTERNS:
        full_pattern = os.path.join(source_dir, pattern)
        for filepath in sorted(glob.glob(full_pattern)):
            rel = os.path.relpath(filepath, source_dir)
            try:
                with open(filepath) as f:
                    content = f.read()
            except UnicodeDecodeError:
                continue
            except IOError as e:
                print(f"Warning: could not read {rel}: {e}", file=sys.stderr)
                continue
            for match in ANNOTATION_RE.finditer(content):
                kv_str = match.group(1)
                fields = dict(KV_RE.findall(kv_str))
                for req in ("observable", "geometry"):
                    if req not in fields:
                        print(f"Warning: COVERAGE annotation in {rel} missing required field '{req}'", file=sys.stderr)
                ann = {
                    "file": rel,
                    "observable": fields.get("observable", ""),
                    "geometry": fields.get("geometry", ""),
                    "material": fields.get("material"),
                    "ref": fields.get("ref"),
                }
                annotations.append(ann)
    return annotations


def cell_key(observable: str, geometry: str, material: str | None) -> tuple[str, str, str]:
    return (observable, geometry, material or "")


def cross_reference(
    cells: list[dict[str, Any]],
    annotations: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    universe_keys: dict[tuple[str, str, str], dict[str, Any]] = {}
    for c in cells:
        k = cell_key(c["observable"], c["geometry"], c["material"])
        if k in universe_keys:
            print(f"Warning: duplicate universe cell {k}", file=sys.stderr)
        universe_keys[k] = c

    coverage: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for ann in annotations:
        k = cell_key(ann["observable"], ann["geometry"], ann["material"])
        if k not in coverage:
            coverage[k] = []
        coverage[k].append(ann)

    results: list[dict[str, Any]] = []
    for c in cells:
        k = cell_key(c["observable"], c["geometry"], c["material"])
        covering = coverage.get(k, [])
        if covering:
            has_ref = any(a.get("ref") for a in covering)
            status = "green" if has_ref else "yellow"
        else:
            status = "red"
        results.append({
            "cell": c,
            "status": status,
            "covering_annotations": covering,
        })

    orphans = []
    for ann in annotations:
        k = cell_key(ann["observable"], ann["geometry"], ann["material"])
        if k not in universe_keys:
            orphans.append(ann)

    return results, orphans


def format_heat_map(results: list[dict[str, Any]]) -> str:
    by_geom: dict[str, list[dict[str, Any]]] = {}
    for r in results:
        geom = r["cell"]["geometry"]
        by_geom.setdefault(geom, []).append(r)

    lines: list[str] = ["## Physics Coverage Heat Map\n"]
    status_icon = {"green": "G", "yellow": "Y", "red": "R"}

    for geom in sorted(by_geom.keys()):
        rows = by_geom.get(geom, [])
        if not rows:
            continue
        lines.append(f"### {geom}\n")
        lines.append("| Observable | Material | Status | Reference |")
        lines.append("|---|---|---|---|")
        for r in sorted(rows, key=lambda x: (x["cell"]["observable"],
                                              x["cell"]["material"] or "")):
            c = r["cell"]
            icon = status_icon.get(r["status"], "?")
            ref = c.get("reference") or (r["covering_annotations"][0].get("ref", "")
                      if r["covering_annotations"] else "")
            mat = c["material"] or "N/A"
            ref_display = ref if ref else "-"
            lines.append(f"| {c['observable']} | {mat} | {icon} | {ref_display} |")
        lines.append("")

    lines.append("**Legend:** G = covered with published reference, "
                 "Y = covered (regression-only), R = uncovered\n")
    return "\n".join(lines)


def format_traceability(results: list[dict[str, Any]]) -> str:
    lines: list[str] = ["## Infrastructure Traceability\n"]
    lines.append("| Observable | Geometry | Material | Test File | Reference |")
    lines.append("|---|---|---|---|---|")

    for r in sorted(results, key=lambda x: (x["cell"]["observable"],
                                             x["cell"]["geometry"],
                                             x["cell"]["material"] or "")):
        c = r["cell"]
        mat = c["material"] or "N/A"
        for ann in r["covering_annotations"]:
            ref = ann.get("ref") or c.get("reference") or "-"
            lines.append(
                f"| {c['observable']} | {c['geometry']} | {mat} "
                f"| `{ann['file']}` | {ref} |")
        if not r["covering_annotations"]:
            lines.append(
                f"| {c['observable']} | {c['geometry']} | {mat} "
                f"| - | - |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    source_dir = os.path.abspath(args.source_dir)

    cells, metadata = load_universe(source_dir)
    annotations = scan_annotations(source_dir)
    results, orphans = cross_reference(cells, annotations)

    n_green = sum(1 for r in results if r["status"] == "green")
    n_yellow = sum(1 for r in results if r["status"] == "yellow")
    n_red = sum(1 for r in results if r["status"] == "red")
    n_orphan = len(orphans)

    red_required = [r for r in results
                    if r["status"] == "red" and r["cell"]["tier"] == "required"]
    red_aspirational = [r for r in results
                        if r["status"] == "red" and r["cell"]["tier"] == "aspirational"]

    report_lines = ["# Validation Coverage Report\n"]
    report_lines.append(f"**Cells:** {len(cells)} total "
                        f"({n_green} green, {n_yellow} yellow, {n_red} red, "
                        f"{n_orphan} orphan)\n")
    report_lines.append(format_heat_map(results))
    report_lines.append(format_traceability(results))

    if orphans:
        report_lines.append("## Orphan Annotations\n")
        report_lines.append("The following annotations reference cells not in the universe file:\n")
        for ann in orphans:
            report_lines.append(
                f"- `{ann['file']}`: observable={ann['observable']} "
                f"geometry={ann['geometry']} material={ann['material'] or 'N/A'}")
        report_lines.append("")

    report = "\n".join(report_lines)

    if args.output:
        with open(args.output, "w") as f:
            f.write(report)
    else:
        print(report)

    summary_lines = [
        f"Coverage: {n_green} green, {n_yellow} yellow, {n_red} red, "
        f"{n_orphan} orphan",
    ]
    if red_required:
        summary_lines.append("UNCOVERED REQUIRED CELLS:")
        for r in red_required:
            c = r["cell"]
            summary_lines.append(
                f"  - {c['observable']}/{c['geometry']}/{c['material'] or 'N/A'}")
    if red_aspirational:
        summary_lines.append(f"Aspirational gaps: {len(red_aspirational)} cells")
    if orphans:
        summary_lines.append(f"Orphan annotations: {n_orphan}")

    for line in summary_lines:
        print(line, file=sys.stderr)

    sys.exit(1 if red_required else 0)


if __name__ == "__main__":
    main()
