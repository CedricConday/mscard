#!/usr/bin/env python3
"""Cohort table and aggregates from an mscard output tree, as Markdown for the README."""

import json
import statistics
import sys
from pathlib import Path


def main(root: Path) -> None:
    rows = []
    for p in sorted(root.iterdir(), key=lambda x: (len(x.name), x.name)):
        m = p / "measure" / "mscard.json"
        g = p / "calibration" / "grade.json"
        if not m.exists():
            continue
        M = json.loads(m.read_text())
        G = json.loads(g.read_text()) if g.exists() else None
        rows.append((p.name, M, G))
    print("| Subject | Years | Baseline lesions | New | Enlarging | Shrinking | Resolved | SEL cand. | Lesion vol %/yr | Brain vol %/yr (eroded) | Boundary | Injected found | False new | Unchanged called changed | SEL/lesion on unchanged | Atrophy recovery |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|")
    agg = {"found": 0, "injected": 0, "false_new": 0, "changed": 0, "unchanged": 0, "resolved": 0, "sel": [], "rec": [], "bound": [], "lights": {}}
    for name, M, G in rows:
        L, br = M["lesions"], M["brain"]["jacobian"]
        gl = (G or {}).get("lights", {})
        n = (G or {}).get("new_lesions", {})
        u = (G or {}).get("unchanged_lesions", {})
        b = (G or {}).get("brain", {})
        s = (G or {}).get("sel", {})
        if G:
            agg["found"] += n["reported"]; agg["injected"] += n["injected"]; agg["false_new"] += n["false_new"]
            agg["changed"] += u["called_enlarging_or_shrinking"]; agg["unchanged"] += u["n"]; agg["resolved"] += u["called_resolved"]
            agg["sel"].append(s["per_baseline_lesion"]); agg["rec"].append(b["jacobian_recovery"])
            for k, v in gl.items():
                agg["lights"].setdefault(k, {}).setdefault(v, 0)
                agg["lights"][k][v] += 1
        agg["bound"].append(br.get("boundary_light"))
        f = lambda x, d=1: "–" if x is None else f"{x:.{d}f}"
        print(f"| {name} | {f(M['interval_years'], 2)} | {L['baseline_count']} | {L['new']} | {L['enlarging']} | {L['shrinking']} | {L['resolved']} | {M['sel']['candidates']} | "
              f"{f(L['volume_change_pct_per_year'])} | {f(br['change_pct_per_year'], 2)} ({f(br.get('change_pct_per_year_eroded'), 2)}) | {br.get('boundary_light', '–')} | "
              f"{n.get('reported', '–')}/{n.get('injected', '–')} | {n.get('false_new', '–')} | {u.get('called_enlarging_or_shrinking', '–')}/{u.get('n', '–')} | "
              f"{f(s.get('per_baseline_lesion'), 2)} | {f(b.get('jacobian_recovery'), 2)} |")
    print()
    bound = {k: agg["bound"].count(k) for k in ("green", "amber", "red")}
    print(f"Subjects: {len(rows)}. Injected new lesions reported: {agg['found']} of {agg['injected']}; false new: {agg['false_new']}. "
          f"Unchanged lesions called enlarging or shrinking: {agg['changed']} of {agg['unchanged']}; called resolved: {agg['resolved']}. "
          f"SEL candidates per unchanged lesion: median {statistics.median(agg['sel']):.2f} (range {min(agg['sel']):.2f} to {max(agg['sel']):.2f}). "
          f"Atrophy recovery: median {statistics.median(agg['rec']):.2f} (range {min(agg['rec']):.2f} to {max(agg['rec']):.2f}). "
          f"Boundary check on the real pairs: {bound}.")
    print("Lights:", json.dumps(agg["lights"]))


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "derivatives/mslesseg"))
