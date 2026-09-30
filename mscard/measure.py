"""Measure one subject: lesiontrack on every follow-up, plus brain volume change from the same warp."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from . import __version__
from .spec import Scan, check_scan_grid, ensure_brainmask, to_lesiontrack


def _vol_mm3(path: Path) -> float:
    img = nib.load(path)
    return float((np.asarray(img.dataobj) > 0).sum() * np.prod(img.header.get_zooms()[:3]))


def brain_volume_change(jacobian: Path, brainmask_halfway: Path | None, baseline_t1_halfway: Path,
                        erode_mm: float = 6.0) -> dict:
    """Whole-brain volume change from the deformation field: mean det J over the baseline brain.

    det J > 1 where the follow-up is locally larger than the baseline, so
    (mean - 1) * 100 is the percent change of brain volume from baseline to
    follow-up, the same quantity lesiontrack integrates inside lesions. The
    same mean over the brain mask eroded by ``erode_mm`` is a self-check the
    calibration cannot provide: two visits stripped differently move the
    boundary, and a change that lives only in the outer rim is a mask
    difference, not tissue.
    """
    from scipy import ndimage as ndi

    img = nib.load(jacobian)
    jac = np.asarray(img.dataobj, dtype=np.float32)
    if brainmask_halfway is not None and Path(brainmask_halfway).exists():
        brain = np.asarray(nib.load(brainmask_halfway).dataobj) > 0
        src = "baseline brain mask"
    else:
        brain = np.asarray(nib.load(baseline_t1_halfway).dataobj) > 0
        src = "nonzero baseline T1w"
    ok = np.isfinite(jac)
    inside = jac[brain & ok]
    depth = ndi.distance_transform_edt(brain, sampling=[float(z) for z in img.header.get_zooms()[:3]])
    core = (depth > erode_mm) & ok
    inner = jac[core]
    full = float((inside.mean() - 1.0) * 100.0)
    er = float((inner.mean() - 1.0) * 100.0) if inner.size else float("nan")
    return {
        "mean_jacobian": float(inside.mean()),
        "change_pct": full,
        "median_jacobian": float(np.median(inside)),
        "brain_voxels": int(brain.sum()),
        "mask_source": src,
        "erode_mm": erode_mm,
        "change_pct_eroded": er,
        "boundary_disagreement_pct": abs(full - er) if np.isfinite(er) else None,
    }


def _fingerprint(scans: list[Scan]) -> dict:
    rec = {}
    for s in scans:
        for what in ("t1", "flair", "mask", "brainmask"):
            q = getattr(s, what)
            if q is not None:
                st = Path(q).stat()
                rec[f"{s.session}/{what}"] = {"path": str(q), "size": st.st_size, "mtime": int(st.st_mtime)}
    return rec


def measure(subject: str, scans: list[Scan], out_dir: Path, threads: int = 4, engine: str = "greedy",
            log: list[str] | None = None, segmenter: str = "given", reg_profile: str = "default") -> dict:
    """Run the longitudinal measurement for one subject and write ``out_dir/mscard.json``."""
    from lesiontrack.config import Params as LTParams
    from lesiontrack.config import RegParams
    from lesiontrack.pipeline import run_subject

    from .figures import overview

    log = log if log is not None else []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if len(scans) < 2:
        raise ValueError(f"{subject}: need a baseline and at least one follow-up")
    scans = sorted(scans, key=lambda s: s.time_years)
    for s in scans:
        check_scan_grid(s)
    scans = [ensure_brainmask(s, out_dir / "brainmask" / f"{s.session}_brainmask.nii.gz", log) for s in scans]
    tps = [to_lesiontrack(s) for s in scans]
    lt_dir = out_dir / "lesiontrack"
    # lesiontrack keys its registration cache on parameters, not inputs: drop it when an input changed.
    fp_path, fp = out_dir / "inputs.json", _fingerprint(scans)
    if fp_path.exists() and json.loads(fp_path.read_text()) != fp and (lt_dir / "reg").exists():
        import shutil

        shutil.rmtree(lt_dir / "reg")
        log.append(f"{subject}: inputs changed since the last run; registration cache dropped")
    fp_path.write_text(json.dumps(fp, indent=2))
    from dataclasses import replace

    reg = RegParams(threads=threads, engine=engine)
    if reg_profile == "smooth":
        # study/sel_sensitivity.py: keeps every injected expansion detected and cuts SEL false positives from 0.50
        # to 0.20 (25 MSLesSeg patients), but under-reads expansion size by about a sixth
        reg = replace(reg, deform_sigma_update="3vox", deform_sigma_total="1vox")
    elif reg_profile != "default":
        raise ValueError(f"unknown registration profile {reg_profile!r}")
    res = run_subject(subject, tps, lt_dir, LTParams(reg=reg))
    for p in res.pairs:
        overview(lt_dir, p.follow_up, out_dir / f"overview_{p.follow_up}.png")
    summary = json.loads((lt_dir / "summary.json").read_text())

    baseline, last = scans[0], scans[-1]
    pair = res.pairs[-1]
    bvc = brain_volume_change(pair.jacobian, pair.baseline_brainmask, pair.baseline_t1)
    dt = pair.dt_years
    bvc["change_pct_per_year"] = bvc["change_pct"] / dt
    bvc["change_pct_per_year_eroded"] = bvc["change_pct_eroded"] / dt
    d = bvc["boundary_disagreement_pct"]
    bvc["boundary_disagreement_pct_per_year"] = None if d is None else d / dt
    dy = bvc["boundary_disagreement_pct_per_year"]
    bvc["boundary_light"] = "none" if dy is None else "green" if dy <= 0.5 else "amber" if dy <= 1.5 else "red"

    per_scan = []
    for s in scans:
        rec = {"session": s.session, "time_years": s.time_years, "date": s.date, "scanner": s.scanner,
               "lesion_volume_mm3": _vol_mm3(s.mask), "lesion_mask": str(s.mask),
               "brain_volume_mm3": _vol_mm3(s.brainmask) if s.brainmask else None}
        per_scan.append(rec)
    mask_bvc = None
    if per_scan[0]["brain_volume_mm3"] and per_scan[-1]["brain_volume_mm3"]:
        r = per_scan[-1]["brain_volume_mm3"] / per_scan[0]["brain_volume_mm3"]
        mask_bvc = {"change_pct": (r - 1) * 100, "change_pct_per_year": (r - 1) * 100 / dt,
                    "note": "ratio of the two brain-mask volumes as given or estimated; independent of the warp, "
                            "but a mask-boundary measure, not a tissue measure"}

    track = res.tracking[res.tracking["follow_up"] == pair.follow_up] if len(res.tracking) else res.tracking
    empty = not len(track) or "class" not in track
    classes = track["class"].value_counts().to_dict() if not empty else {}
    lv0 = float(track["volume_baseline_mm3"].sum()) if not empty else 0.0
    lv1 = float(track["volume_followup_mm3"].sum()) if not empty else 0.0
    cand = res.candidates
    rows = track.sort_values("volume_followup_mm3", ascending=False).to_dict("records") if not empty else []

    rec = {
        "mscard_version": __version__,
        "lesiontrack_version": summary["lesiontrack_version"],
        "subject": subject,
        "segmenter": segmenter,
        "reg_profile": reg_profile,
        "baseline": baseline.session,
        "follow_up": last.session,
        "follow_ups": summary["follow_ups"],
        "interval_years": dt,
        "voxel_mm3": summary["voxel_mm3"],
        "scans": per_scan,
        "lesions": {
            "baseline_count": summary["baseline_lesion_count"],
            "baseline_volume_mm3": summary["baseline_lesion_volume_mm3"],
            "followup_volume_mm3": lv1,
            "volume_change_pct": (lv1 - lv0) / lv0 * 100 if lv0 else None,
            "volume_change_pct_per_year": (lv1 - lv0) / lv0 * 100 / dt if lv0 else None,
            "new": int(classes.get("new", 0)),
            "adjacent_fragment": int(classes.get("adjacent_fragment", 0)),
            "enlarging": int(classes.get("enlarging", 0)),
            "shrinking": int(classes.get("shrinking", 0)),
            "resolved": int(classes.get("resolved", 0)),
            "stable": int(classes.get("stable", 0)),
            "trend_up": int(classes.get("trend_up", 0)),
            "trend_down": int(classes.get("trend_down", 0)),
            "groups": len(rows),
            "table": rows,
        },
        "sel": {
            "candidates": len(cand),
            "candidate_volume_mm3": float(cand["volume_mm3"].sum()) if len(cand) else 0.0,
            "definite": None,
            "note": "definite SELs are a cohort-relative ranking (z-scores across subjects, Elliott 2019); "
                    "run `lesiontrack cohort` over several subjects to obtain it",
            "table": cand.to_dict("records") if len(cand) else [],
        },
        "brain": {"jacobian": bvc, "mask_ratio": mask_bvc},
        "images": {"overview": str(out_dir / f"overview_{pair.follow_up}.png")},
        "params": summary["params"],
        "log": log,
    }
    from .calibrate import scrub

    rec = scrub(rec)
    with open(out_dir / "mscard.json", "w") as fh:
        json.dump(rec, fh, indent=2, allow_nan=False)
    return rec

