"""Calibration: a synthetic follow-up with known truth, the same pipeline, and a grade.

The synthetic follow-up is the subject's own baseline scan, one modelled year
later: the whole brain contracted by a known factor (bidsgate inject-atrophy),
a known set of new lesions placed away from the existing ones (bidsgate
inject-lesions), a small rigid repositioning, a gain change and noise. Every
baseline lesion is otherwise unchanged, so the truth for each report line is
known: which follow-up lesions are new, that no baseline lesion enlarged,
shrank, resolved or expanded slowly, and the brain volume change.
"""

from __future__ import annotations

import json
import zlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .config import Bands, SyntheticParams
from .spec import Scan

SESSION = "SYN"


def _save(arr: np.ndarray, like: nib.Nifti1Image, path: Path, dtype=None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(arr.astype(dtype) if dtype is not None else arr, like.affine), path)
    return path


def make_synthetic(baseline: Scan, out_dir: Path, sp: SyntheticParams, log: list[str] | None = None) -> tuple[Scan, dict]:
    """Write the synthetic follow-up of ``baseline`` under ``out_dir`` and return (scan, truth)."""
    from bidsgate.inject_atrophy import _warp, displacement
    from bidsgate.inject_lesions import LesionSpec, inject
    from lesiontrack.synth import _rigid_field

    out_dir = Path(out_dir)
    truth_json = out_dir / "truth.json"
    t1_p, fl_p = out_dir / "T1w.nii.gz", out_dir / "FLAIR.nii.gz"
    mask_p, brain_p, new_p = out_dir / "mask.nii.gz", out_dir / "brainmask.nii.gz", out_dir / "newlesions_truth.nii.gz"
    if truth_json.exists() and all(p.exists() for p in (t1_p, fl_p, mask_p, brain_p, new_p)):
        truth = json.loads(truth_json.read_text())
        if truth.get("params") == sp.__dict__ and truth.get("source", {}).get("mask") == str(baseline.mask):
            return _scan(baseline, out_dir, sp), truth
    if baseline.mask is None or baseline.brainmask is None:
        raise ValueError("make_synthetic needs a baseline with a lesion mask and a brain mask")
    seed = (sp.seed * 1_000_003 + zlib.crc32(baseline.subject.encode())) % (2**32)
    rng = np.random.default_rng(seed)
    img = nib.load(baseline.t1)
    zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
    t1 = np.asarray(img.dataobj, dtype=np.float32)
    fl = np.asarray(nib.load(baseline.flair).dataobj, dtype=np.float32)
    mask = np.asarray(nib.load(baseline.mask).dataobj) > 0
    brain = np.asarray(nib.load(baseline.brainmask).dataobj) > 0

    # 1. Whole-brain contraction by a known factor; every image and mask follows the same field.
    disp = displacement(t1.shape, zooms, brain, sp.volume_factor, sp.falloff_mm)
    t1a = _warp(t1, disp, 1)
    fla = _warp(fl, disp, 1)
    maska = _warp(mask.astype(np.uint8), disp, 0) > 0
    braina = _warp(brain.astype(np.uint8), disp, 0) > 0
    del disp

    # 2. New lesions, placed by bidsgate inside a brain mask that excludes a margin around every
    #    baseline lesion, so each injected lesion is truly new rather than an enlargement.
    tmp = out_dir / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    it = round(sp.exclusion_mm / min(zooms))
    excl = ndi.binary_dilation(maska, structure=ndi.generate_binary_structure(3, 2), iterations=max(1, it))
    place = _save((braina & ~excl).astype(np.uint8), img, tmp / "placement.nii.gz")
    _save(t1a, img, tmp / "T1w_atrophied.nii.gz")
    _save(fla, img, tmp / "FLAIR_atrophied.nii.gz")
    new_json = out_dir / "newlesions_truth.json"
    spec = LesionSpec(n=sp.n_new, seed=int(rng.integers(0, 2**31 - 1)))
    inject(tmp / "T1w_atrophied.nii.gz", tmp / "FLAIR_atrophied.nii.gz", tmp / "T1w_injected.nii.gz",
           tmp / "FLAIR_injected.nii.gz", tmp / "newlesions.nii.gz", new_json, spec, mask_path=place)
    t1b = np.asarray(nib.load(tmp / "T1w_injected.nii.gz").dataobj, dtype=np.float32)
    flb = np.asarray(nib.load(tmp / "FLAIR_injected.nii.gz").dataobj, dtype=np.float32)
    newb = np.asarray(nib.load(tmp / "newlesions.nii.gz").dataobj).astype(np.int32)
    new_truth = json.loads(new_json.read_text())

    # 3. Rigid repositioning, gain, noise: what a second visit adds even when nothing changed.
    rot = rng.uniform(-1, 1, 3) * sp.rigid_rot_deg
    trans = rng.uniform(-1, 1, 3) * sp.rigid_trans_vox
    rig = _rigid_field(t1.shape, rot, trans)
    t1c = _warp(t1b, rig, 1)
    flc = _warp(flb, rig, 1)
    maskc = _warp((maska | (newb > 0)).astype(np.uint8), rig, 0) > 0
    newc = _warp(newb, rig, 0)
    brainc = _warp(braina.astype(np.uint8), rig, 0) > 0
    del rig
    gains = rng.uniform(*sp.intensity_scale, 2)
    head = t1c > 0
    for arr, g in ((t1c, gains[0]), (flc, gains[1])):
        arr *= g
        sd = sp.noise_frac * float(np.percentile(arr[head], 99)) if head.any() else 0.0
        arr[head] += rng.normal(0, sd, int(head.sum())).astype(np.float32)
        np.maximum(arr, 0, out=arr)
    _save(t1c, img, t1_p)
    _save(flc, img, fl_p)
    _save(maskc, img, mask_p, np.uint8)
    _save(brainc, img, brain_p, np.uint8)
    _save(newc, img, new_p, np.int32)
    for p in tmp.glob("*"):
        p.unlink()
    tmp.rmdir()
    vox = float(np.prod(zooms))
    truth = {
        "kind": "mscard synthetic follow-up",
        "session": SESSION,
        "dt_years": sp.dt_years,
        "volume_factor": sp.volume_factor,
        "injected_brain_change_pct": (sp.volume_factor - 1) * 100,
        "n_new": sp.n_new,
        "new_lesions": new_truth["lesions"],
        "new_lesion_volume_mm3": new_truth["total_volume_mm3"],
        "baseline_lesion_count_truth": int(ndi.label(maska, structure=ndi.generate_binary_structure(3, 2))[1]),
        "baseline_lesion_volume_after_contraction_mm3": float(maska.sum() * vox),
        "rigid": {"rot_deg": rot.tolist(), "trans_vox": trans.tolist()},
        "gains": gains.tolist(),
        "seed": seed,
        "params": sp.__dict__,
        "source": {"t1": str(baseline.t1), "flair": str(baseline.flair), "mask": str(baseline.mask),
                   "brainmask": str(baseline.brainmask)},
        "truth_per_line": {
            "new": "exactly the injected lesions; any other follow-up-only component is false",
            "enlarging/shrinking/resolved": "none; every baseline lesion is unchanged apart from the global factor",
            "sel": "none; no lesion expands",
            "brain_volume_change_pct": (sp.volume_factor - 1) * 100,
        },
    }
    truth_json.write_text(json.dumps(truth, indent=2))
    if log is not None:
        log.append(f"{baseline.subject}: synthetic follow-up, factor {sp.volume_factor}, {sp.n_new} new lesions, seed {seed}")
    return _scan(baseline, out_dir, sp), truth


def _scan(baseline: Scan, out_dir: Path, sp: SyntheticParams) -> Scan:
    return Scan(subject=baseline.subject, session=SESSION, time_years=baseline.time_years + sp.dt_years,
                t1=out_dir / "T1w.nii.gz", flair=out_dir / "FLAIR.nii.gz", mask=out_dir / "mask.nii.gz",
                brainmask=out_dir / "brainmask.nii.gz", date="synthetic", scanner=baseline.scanner)


def _ids(v) -> list[int]:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return []
    # pandas reads an all-single-id column as float ("1.0"); a multi-id cell is "1,3"
    return [int(float(x)) for x in str(v).replace(";", ",").split(",") if str(x).strip()]


def light(value, green, amber, kind: str) -> str:
    """kind: 'min' (higher is better), 'max' (lower is better), 'range' (inside is better)."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "none"
    if kind == "min":
        return "green" if value >= green else "amber" if value >= amber else "red"
    if kind == "max":
        return "green" if value <= green else "amber" if value <= amber else "red"
    lo, hi = green
    alo, ahi = amber
    return "green" if lo <= value <= hi else "amber" if alo <= value <= ahi else "red"


def grade(baseline: Scan, syn: Scan, cal_dir: Path, measurement: dict, truth: dict, segmenter: str,
          bands: Bands, threads: int = 4, reference_mask: Path | None = None) -> dict:
    """Compare the calibration measurement with the injected truth, line by line."""
    from lesiontrack.registration import _reslice

    cal_dir = Path(cal_dir)
    lt = cal_dir / "measure" / "lesiontrack"
    fu = measurement["follow_up"]
    half = lt / "reg" / f"{fu}_halfway.mat"
    truth_half = cal_dir / "newlesions_truth_halfway.nii.gz"
    _reslice(baseline.t1, syn.t1.parent / "newlesions_truth.nii.gz", truth_half, f"{half}", "NN", threads,
             cal_dir / "reslice.log")
    tn = np.asarray(nib.load(truth_half).dataobj).astype(np.int32)
    fl = np.asarray(nib.load(lt / f"lesion_labels_{fu}.nii.gz").dataobj).astype(np.int32)
    table = pd.read_csv(lt / "lesion_tracking.tsv", sep="\t", dtype={"baseline_ids": str, "followup_ids": str, "follow_up": str})
    table = table[table["follow_up"].astype(str) == fu]
    new_ids = [i for _, r in table[table["class"] == "new"].iterrows() for i in _ids(r["followup_ids"])]
    frag_ids = [i for _, r in table[table["class"] == "adjacent_fragment"].iterrows() for i in _ids(r["followup_ids"])]
    pred_new = np.isin(fl, new_ids)
    pred_frag = np.isin(fl, frag_ids)
    seg = fl > 0
    per = []
    for les in truth["new_lesions"]:
        sel = tn == les["id"]
        n = int(sel.sum())
        per.append({
            "id": les["id"], "volume_mm3": les["volume_mm3"], "voxels_halfway": n,
            "segmented": bool(seg[sel].any()) if n else False,
            "reported_new": bool(pred_new[sel].any()) if n else False,
            "reported_fragment": bool(pred_frag[sel].any()) if n else False,
        })
    dfp = pd.DataFrame(per)
    dfp["reported"] = dfp["reported_new"] | dfp["reported_fragment"]
    tn_near = ndi.binary_dilation(tn > 0, structure=ndi.generate_binary_structure(3, 2), iterations=1)
    false_new = [i for i in new_ids if not tn_near[fl == i].any()]
    false_frag = [i for i in frag_ids if not tn_near[fl == i].any()]
    by_size = {}
    for lo, hi, name in ((0, 50, "<50"), (50, 150, "50-150"), (150, 400, "150-400"), (400, 1e9, ">=400")):
        s = dfp[(dfp["volume_mm3"] >= lo) & (dfp["volume_mm3"] < hi)]
        if len(s):
            by_size[name] = {"n": len(s), "reported": int(s["reported"].sum()), "segmented": int(s["segmented"].sum())}

    base = table[table["n_baseline"] > 0]
    n_base = len(base)
    cls = base["class"].value_counts().to_dict()
    false_change = int(cls.get("enlarging", 0) + cls.get("shrinking", 0))
    false_trend = int(cls.get("trend_up", 0) + cls.get("trend_down", 0))
    false_resolved = int(cls.get("resolved", 0))

    inj = truth["injected_brain_change_pct"]
    jac = measurement["brain"]["jacobian"]["change_pct"]
    mr = measurement["brain"]["mask_ratio"]["change_pct"] if measurement["brain"].get("mask_ratio") else None
    n_cand = measurement["sel"]["candidates"]
    base_count = measurement["lesions"]["baseline_count"] or 1

    seg_grade = None
    if segmenter == "lst-ai":
        seg_grade = _segmentation_grade(baseline, syn, cal_dir, truth, reference_mask)

    b = bands
    g = {
        "segmenter": segmenter,
        "truth": {k: truth[k] for k in ("volume_factor", "injected_brain_change_pct", "n_new", "new_lesion_volume_mm3", "rigid", "gains", "seed")},
        "new_lesions": {
            "injected": len(dfp),
            "reported_as_new": int(dfp["reported_new"].sum()),
            "reported_as_fragment": int(dfp["reported_fragment"].sum()),
            "reported": int(dfp["reported"].sum()),
            "segmented": int(dfp["segmented"].sum()),
            "sensitivity": float(dfp["reported"].mean()) if len(dfp) else None,
            "smallest_reported_mm3": float(dfp[dfp["reported"]]["volume_mm3"].min()) if dfp["reported"].any() else None,
            "largest_missed_mm3": float(dfp[~dfp["reported"]]["volume_mm3"].max()) if (~dfp["reported"]).any() else None,
            "false_new": len(false_new),
            "false_fragment": len(false_frag),
            "by_size": by_size,
            "per_lesion": per,
        },
        "unchanged_lesions": {
            "n": n_base,
            "called_enlarging_or_shrinking": false_change,
            "called_trend": false_trend,
            "called_resolved": false_resolved,
            "false_change_rate": false_change / n_base if n_base else None,
            "false_resolved_rate": false_resolved / n_base if n_base else None,
            "classes": {k: int(v) for k, v in cls.items()},
        },
        "sel": {"candidates": n_cand, "per_baseline_lesion": n_cand / base_count, "truth": 0},
        "brain": {
            "injected_change_pct": inj,
            "jacobian_change_pct": jac,
            "jacobian_recovery": jac / inj if inj else None,
            "mask_ratio_change_pct": mr,
            "mask_ratio_recovery": (mr / inj) if (mr is not None and inj) else None,
        },
        "segmentation": seg_grade,
    }
    g["lights"] = {
        "new_sensitivity": light(g["new_lesions"]["sensitivity"], *b.new_sensitivity_min, "min"),
        "false_new": light(g["new_lesions"]["false_new"], *b.false_new_max, "max"),
        "atrophy_recovery": light(g["brain"]["jacobian_recovery"], *b.atrophy_recovery_range, "range"),
        "false_change_rate": light(g["unchanged_lesions"]["false_change_rate"], *b.false_change_rate_max, "max"),
        "false_resolved_rate": light(g["unchanged_lesions"]["false_resolved_rate"], *b.false_resolved_rate_max, "max"),
        "sel_candidates_per_lesion": light(g["sel"]["per_baseline_lesion"], *b.sel_candidates_per_lesion_max, "max"),
        "seg_sensitivity": light(seg_grade["sensitivity"] if seg_grade else None, *b.seg_sensitivity_min, "min"),
        "seg_false_positives": light(seg_grade["false_positive_components"] if seg_grade else None, *b.seg_false_positive_max, "max"),
    }
    g["bands"] = bands.__dict__
    with open(cal_dir / "grade.json", "w") as fh:
        json.dump(g, fh, indent=2, default=_default)
    return g


def _segmentation_grade(baseline: Scan, syn: Scan, cal_dir: Path, truth: dict, reference_mask: Path | None) -> dict:
    """LST-AI on the synthetic follow-up against the injected lesions, old lesions masked out."""
    from bidsgate.score import score_lesions

    new_truth = syn.t1.parent / "newlesions_truth.nii.gz"
    old = np.asarray(nib.load(syn.t1.parent / "mask.nii.gz").dataobj) > 0
    newl = np.asarray(nib.load(new_truth).dataobj) > 0
    old_only = old & ~newl
    img = nib.load(syn.mask)
    pred = np.asarray(img.dataobj) > 0
    zooms = img.header.get_zooms()[:3]
    near_old = ndi.distance_transform_edt(~old_only, sampling=zooms) <= 2.0
    pred_new_only = pred & ~near_old
    p = cal_dir / "lstai_pred_newonly.nii.gz"
    nib.save(nib.Nifti1Image(pred_new_only.astype(np.uint8), img.affine), p)
    sc = score_lesions(syn.t1.parent / "newlesions_truth.json", new_truth, p)
    out = {k: sc[k] for k in ("dice", "sensitivity", "lesions", "detected", "false_positive_components",
                              "false_positive_volume_mm3", "by_size", "volume_ratio")}
    out["note"] = "LST-AI mask of the synthetic follow-up scored against the injected lesions only; " \
                  "predicted voxels within 2 mm of a baseline lesion are excluded first"
    if reference_mask is not None and baseline.mask is not None:
        a = np.asarray(nib.load(reference_mask).dataobj) > 0
        b = np.asarray(nib.load(baseline.mask).dataobj) > 0
        out["baseline_reference_dice"] = float(2 * (a & b).sum() / (a.sum() + b.sum())) if (a.sum() + b.sum()) else None
        out["baseline_reference_mask"] = str(reference_mask)
    return out


def _default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(str(type(o)))
