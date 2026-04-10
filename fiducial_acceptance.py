#!/usr/bin/env python3
#python3 fiducial_acceptance.py
"""
Compute fiducial acceptances from files in fidXS/.
(A bit porcata, but that's what the convento is passing)

Formula used:
    fidXS = XS * BR * Acc * 1000
    -> Acc = fidXS / (XS * BR * 1000)
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path


# Keep these easy to edit.
XS = {
    "ggH": 51.96,
    "VBFH": 4.067,
    "VH": 2.3781,
    "ttH": 0.5638,
}
BR = 0.00227
OUTPUT_TXT = "fiducial_acceptance.txt"

TARGET_OBS = {"NJ", "PTH", "PTJ0"}
TARGET_MODES = {"ggH", "VBFH", "VH", "ttH"}
FILE_RE = re.compile(r"^fidXS_(?P<obs>NJ|PTH|PTJ0)_(?P<mode>ggH|VBFH|VH|ttH)\.py$")


def load_python_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load file: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compute_acceptance(fid_xs_values, xs_pb: float, br: float):
    denom = xs_pb * br * 1000.0
    return [value / denom for value in fid_xs_values]


xs_map = dict(XS)
br = BR
folder = Path("fidXS")
if not folder.exists():
    raise FileNotFoundError(f"Folder not found: {folder}")

results = {}
files = sorted(folder.glob("fidXS_*.py"))
for path in files:
    match = FILE_RE.match(path.name)
    if not match:
        continue

    obs = match.group("obs")
    mode = match.group("mode")
    if obs not in TARGET_OBS or mode not in TARGET_MODES:
        continue

    module = load_python_module(path)
    if not hasattr(module, "fidXS"):
        print(f"Skipping {path.name}: no fidXS vector")
        continue

    fid_xs = list(module.fidXS)
    acc = compute_acceptance(fid_xs, xs_map[mode], br)
    key = f"{obs}_{mode}"
    results[key] = {
        "file": str(path),
        "xs_pb": xs_map[mode],
        "br": br,
        "fidXS": fid_xs,
        "Acc": acc,
    }

if not results:
    print("No matching files found.")
else:
    lines = []
    for key in sorted(results):
        lines.append(f"[{key}]")
        lines.append(f"Acc = {results[key]['Acc']}")
        lines.append("")

    output_text = "\n".join(lines).rstrip() + "\n"
    print(output_text, end="")
    Path(OUTPUT_TXT).write_text(output_text)
    print(f"\nSaved output to {OUTPUT_TXT}")
