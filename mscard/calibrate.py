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
import shutil
import zlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .config import Bands, SyntheticParams
from .spec import Scan


def label_lesions(mask: np.ndarray, connectivity: int = 2) -> np.ndarray:
    from lesiontrack.tracking import label_lesions as _ll

    return _ll(mask, connectivity)

SESSION = "SYN"
TRUTH_FORMAT = 2  # bump when the synthetic construction or the truth record changes; old caches are rebuilt


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
    params_json = json.loads(json.dumps(sp.__dict__))  # tuples become lists, as in the stored truth
    if truth_json.exists() and all(p.exists() for p in (t1_p, fl_p, mask_p, brain_p, new_p)):
        truth = json.loads(truth_json.read_text())
        if truth.get("format") == TRUTH_FORMAT and truth.get("params") == params_json and truth.get("source") == _source(baseline):
            truth["rewritten"] = False
            return _scan(baseline, out_dir, sp), truth
    if baseline.mask is None or baseline.brainmask is None:
        raise ValueError("make_synthetic needs a baseline with a lesion mask and a brain mask")
    # A stale synthetic invalidates every downstream cache (lesiontrack keys its registration
    # cache on parameters only, not on the input files).
    for stale in (out_dir.parent / "measure", out_dir.parent / "lst-ai", out_dir.parent / "grade.json",
                  out_dir.parent / "newlesions_truth_halfway.nii.gz"):
        if stale.is_dir():
            shutil.rmtree(stale)
        elif stale.exists():
            stale.unlink()
    seed = (sp.seed * 1_000_003 + zlib.crc32(baseline.subject.encode())) % (2**32)
    rng = np.random.default_rng(seed)
    img = nib.load(baseline.t1)
    zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
    if max(zooms) / min(zooms) > 1.05:
        raise ValueError(f"make_synthetic needs isotropic voxels (got {zooms}); the rigid perturbation is built in "
                         "index space and would add a strain of the size of the injected atrophy on anisotropic data. "
                         "Resample to isotropic first.")
    t1 = np.asarray(img.dataobj, dtype=np.float32)
    fl = np.asarray(nib.load(baseline.flair).dataobj, dtype=np.float32)
    mask = np.asarray(nib.load(baseline.mask).dataobj) > 0
    brain = np.asarray(nib.load(baseline.brainmask).dataobj) > 0

    # 1. Whole-brain contraction by a known factor; every image and mask follows the same field.
    disp = displacement(t1.shape, zooms, brain, sp.volume_factor, sp.falloff_mm)
    t1a = _warp(t1, disp, 1)
    fla = _warp(fl, disp, 1)
    # Masks follow the field as the zero level set of their signed distance: nearest-neighbour
    # keeps the old boundary, and a 0.5 threshold on an interpolated binary mask erodes small
    # lesions; the signed distance moves the boundary by the field and nothing else.
    lab_bl = label_lesions(mask, 2)
    laba = _warp_labels(lab_bl, disp, zooms)
    maska = laba > 0
    braina = _warp_mask(brain, disp, zooms)
    del disp

    # 2. New lesions, placed by bidsgate inside a brain mask that excludes a margin around every
    #    baseline lesion, so each injected lesion is truly new rather than an enlargement.
    tmp = out_dir / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    _save(t1a, img, tmp / "T1w_atrophied.nii.gz")
    _save(fla, img, tmp / "FLAIR_atrophied.nii.gz")
    new_json = out_dir / "newlesions_truth.json"
    spec = LesionSpec(n=sp.n_new, seed=int(rng.integers(0, 2**31 - 1)))
    exclusion_used = None
    # A heavily lesioned brain may leave too little white matter once a margin around every lesion is
    # excluded; step the margin down rather than give up, and record the margin that was used.
    for margin in sorted({sp.exclusion_mm, 4.0, 2.0}, reverse=True):
        if margin > sp.exclusion_mm:
            continue
        it = max(1, round(margin / min(zooms)))
        excl = ndi.binary_dilation(maska, structure=ndi.generate_binary_structure(3, 2), iterations=it)
        _save((braina & ~excl).astype(np.uint8), img, tmp / "placement.nii.gz")
        try:
            inject(tmp / "T1w_atrophied.nii.gz", tmp / "FLAIR_atrophied.nii.gz", tmp / "T1w_injected.nii.gz",
                   tmp / "FLAIR_injected.nii.gz", tmp / "newlesions.nii.gz", new_json, spec, mask_path=tmp / "placement.nii.gz")
            exclusion_used = margin
            break
        except ValueError as e:
            if log is not None:
                log.append(f"{baseline.subject}: lesion placement failed with a {margin:g} mm margin ({e}); trying smaller")
    if exclusion_used is None:
        raise ValueError(f"{baseline.subject}: no room to place new lesions even with a 2 mm margin")
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
    labc = _warp_labels(laba, rig, zooms)
    newc = _warp_labels(newb, rig, zooms)
    maskc = (labc > 0) | (newc > 0)
    brainc = _warp_mask(braina, rig, zooms)
    del rig
    vox = float(np.prod(zooms))
    vb, va = np.bincount(lab_bl.ravel()), np.bincount(labc.ravel(), minlength=lab_bl.max() + 1)
    baseline_lesions = [{"id": int(i), "volume_before_mm3": float(vb[i] * vox), "volume_after_mm3": float(va[i] * vox)}
                        for i in range(1, lab_bl.max() + 1)]
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
    _save(labc, img, out_dir / "baselinelesions_truth.nii.gz", np.int32)
    for p in tmp.glob("*"):
        p.unlink()
    tmp.rmdir()
    truth = {
        "rewritten": True,
        "format": TRUTH_FORMAT,
        "kind": "mscard synthetic follow-up",
        "session": SESSION,
        "dt_years": sp.dt_years,
        "volume_factor": sp.volume_factor,
        "injected_brain_change_pct": (sp.volume_factor - 1) * 100,
        "n_new": sp.n_new,
        "exclusion_mm_used": exclusion_used,
        "new_lesions": new_truth["lesions"],
        "new_lesion_volume_mm3": new_truth["total_volume_mm3"],
        "baseline_lesion_count_truth": int(lab_bl.max()),
        "baseline_lesion_volume_mm3": float(mask.sum() * vox),
        "baseline_lesion_volume_after_mm3": float((labc > 0).sum() * vox),
        "baseline_lesions": baseline_lesions,
        "rigid": {"rot_deg": rot.tolist(), "trans_vox": trans.tolist()},
        "gains": gains.tolist(),
        "seed": seed,
        "params": sp.__dict__,
        "source": _source(baseline),
        "truth_per_line": {
            "new": "exactly the injected lesions; any other follow-up-only component is false",
            "enlarging/shrinking/resolved": "per lesion in baseline_lesions: volume_after_mm3 against volume_before_mm3 "
                                           "(the global factor plus voxel quantisation of the moved boundary); no lesion resolves",
            "sel": "none; no lesion expands",
            "brain_volume_change_pct": (sp.volume_factor - 1) * 100,
        },
    }
    truth_json.write_text(json.dumps(scrub(truth), indent=2, allow_nan=False))
    if log is not None:
        log.append(f"{baseline.subject}: synthetic follow-up, factor {sp.volume_factor}, {sp.n_new} new lesions, seed {seed}")
    return _scan(baseline, out_dir, sp), truth


def _source(baseline: Scan) -> dict:
    """The inputs the synthetic was built from, with size and mtime so a changed mask is noticed."""
    rec = {}
    for what in ("t1", "flair", "mask", "brainmask"):
        q = getattr(baseline, what)
        st = Path(q).stat() if q is not None else None
        rec[what] = None if q is None else {"path": str(q), "size": st.st_size, "mtime": int(st.st_mtime)}
    return rec


def _warp_mask(mask: np.ndarray, disp: np.ndarray, zooms: tuple) -> np.ndarray:
    """Move a binary mask by a sampling field as the zero level set of its signed distance (mm)."""
    from bidsgate.inject_atrophy import _warp

    sd = (ndi.distance_transform_edt(~mask, sampling=zooms) - ndi.distance_transform_edt(mask, sampling=zooms)).astype(np.float32)
    return _warp(sd, disp, 1, cval=1e3) <= 0


def _warp_labels(labels: np.ndarray, disp: np.ndarray, zooms: tuple) -> np.ndarray:
    """Warp an integer label map: the support as a level set, labels by nearest neighbour then nearest fill."""
    from bidsgate.inject_atrophy import _warp

    support = _warp_mask(labels > 0, disp, zooms)
    nn = _warp(labels, disp, 0)
    missing = support & (nn == 0)
    if missing.any():
        idx = ndi.distance_transform_edt(nn == 0, return_distances=False, return_indices=True)
        nn = np.where(missing, nn[tuple(idx)], nn)
    return np.where(support, nn, 0).astype(np.int32)


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
    from lesiontrack.registration import run_greedy

    cal_dir = Path(cal_dir)
    lt = cal_dir / "measure" / "lesiontrack"
    fu = measurement["follow_up"]
    if len(measurement["follow_ups"]) != 1:
        raise ValueError("grade expects the calibration measurement to hold exactly one follow-up")
    # Single pair: lesiontrack resliced the follow-up mask with its own halfway matrix (frame=None).
    half = lt / "reg" / f"{fu}_halfway.mat"
    truth_half = cal_dir / "newlesions_truth_halfway.nii.gz"
    run_greedy(f"-d 3 -threads {threads} -rf {baseline.t1} -ri NN -rt short -rm {syn.t1.parent / 'newlesions_truth.nii.gz'} "
               f"{truth_half} -r {half}", cal_dir / "reslice.log")
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
    dfp = pd.DataFrame(per, columns=["id", "volume_mm3", "voxels_halfway", "segmented", "reported_new", "reported_fragment"])
    dfp["reported"] = dfp["reported_new"].astype(bool) | dfp["reported_fragment"].astype(bool)
    tn_near = ndi.binary_dilation(tn > 0, structure=ndi.generate_binary_structure(3, 2), iterations=1)
    false_new = [i for i in new_ids if not tn_near[fl == i].any()]
    false_frag = [i for i in frag_ids if not tn_near[fl == i].any()]
    by_size = {}
    for lo, hi, name in ((0, 50, "<50"), (50, 150, "50-150"), (150, 400, "150-400"), (400, 1e9, ">=400")):
        s = dfp[(dfp["volume_mm3"] >= lo) & (dfp["volume_mm3"] < hi)]
        if len(s):
            by_size[name] = {"n": len(s), "reported": int(s["reported"].sum()), "segmented": int(s["segmented"].sum())}

    # Baseline lesions: the truth for each is its own volume in the synthetic mask (the global
    # factor plus the voxel quantisation of a moved boundary), so a call is judged against the
    # class that truth rate falls in, with lesiontrack's own bands.
    from lesiontrack.config import TrackParams
    from lesiontrack.tracking import classify_rate

    bl_truth_half = cal_dir / "baselinelesions_truth_halfway.nii.gz"
    run_greedy(f"-d 3 -threads {threads} -rf {baseline.t1} -ri NN -rt short -rm {syn.t1.parent / 'baselinelesions_truth.nii.gz'} "
               f"{bl_truth_half} -r {half}", cal_dir / "reslice.log")
    tb = np.asarray(nib.load(bl_truth_half).dataobj).astype(np.int32)
    lb = np.asarray(nib.load(lt / "baseline_lesion_labels.nii.gz").dataobj).astype(np.int32)
    tinfo = {x["id"]: x for x in truth.get("baseline_lesions", [])}
    dt = truth["dt_years"]
    vox = measurement["voxel_mm3"]
    tp = TrackParams()
    base = table[table["n_baseline"] > 0]
    judged = []
    for _, r in base.iterrows():
        sel = np.isin(lb, _ids(r["baseline_ids"]))
        ids = [int(i) for i in np.unique(tb[sel]) if i > 0 and int(i) in tinfo]
        if not ids:
            continue
        before = sum(tinfo[i]["volume_before_mm3"] for i in ids)
        after = sum(tinfo[i]["volume_after_mm3"] for i in ids)
        if after <= 0:
            t_cls, t_rate = "resolved", -100.0 / dt
        else:
            t_rate = (after - before) / before * 100.0 / dt
            t_cls = classify_rate(t_rate, round((after - before) / vox), tp)
        judged.append({"group_id": int(r["group_id"]), "predicted": str(r["class"]), "truth": t_cls,
                       "truth_rate_pct_per_year": t_rate, "predicted_rate_pct_per_year": None if r["class"] == "resolved" else float(r["pct_per_year"]),
                       "n_lesions": len(ids)})
    n_base = sum(j["n_lesions"] for j in judged)
    cls = {}
    for j in judged:
        cls[j["predicted"]] = cls.get(j["predicted"], 0) + 1
    changed = {"enlarging", "shrinking"}
    false_change = sum(1 for j in judged if j["predicted"] in changed and j["truth"] not in changed)
    missed_change = sum(1 for j in judged if j["truth"] in changed and j["predicted"] not in changed and j["predicted"] != "resolved")
    false_trend = sum(1 for j in judged if j["predicted"] in ("trend_up", "trend_down") and j["truth"] == "stable")
    false_resolved = sum(1 for j in judged if j["predicted"] == "resolved" and j["truth"] != "resolved")
    agree = sum(1 for j in judged if j["predicted"] == j["truth"])
    truth_changed = sum(1 for j in judged if j["truth"] in changed)

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
            "groups": len(judged),
            "called_enlarging_or_shrinking": false_change,
            "missed_enlarging_or_shrinking": missed_change,
            "truth_enlarging_or_shrinking": truth_changed,
            "called_trend": false_trend,
            "called_resolved": false_resolved,
            "class_agreement": agree / len(judged) if judged else None,
            "false_change_rate": false_change / len(judged) if judged else None,
            "false_resolved_rate": false_resolved / len(judged) if judged else None,
            "classes": {k: int(v) for k, v in cls.items()},
            "per_group": judged,
            "note": "truth per baseline lesion is its own volume in the synthetic mask; a call is false when it lands "
                    "in the enlarging/shrinking band while the truth rate does not, and resolved is false when the lesion is still there",
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
    # With supplied masks the synthetic follow-up mask contains the injected lesions by construction,
    # so their detection is the annotator's, not the pipeline's: those two lights are not awarded.
    detection_graded = segmenter != "given"
    g["lights"] = {
        "new_sensitivity": light(g["new_lesions"]["sensitivity"], *b.new_sensitivity_min, "min") if detection_graded else "none",
        "false_new": light(g["new_lesions"]["false_new"], *b.false_new_max, "max") if detection_graded else "none",
        "atrophy_recovery": light(g["brain"]["jacobian_recovery"], *b.atrophy_recovery_range, "range"),
        "false_change_rate": light(g["unchanged_lesions"]["false_change_rate"], *b.false_change_rate_max, "max"),
        "false_resolved_rate": light(g["unchanged_lesions"]["false_resolved_rate"], *b.false_resolved_rate_max, "max"),
        "sel_candidates_per_lesion": light(g["sel"]["per_baseline_lesion"], *b.sel_candidates_per_lesion_max, "max"),
        "seg_sensitivity": light(seg_grade["sensitivity"] if seg_grade else None, *b.seg_sensitivity_min, "min"),
        "seg_false_positives": light(seg_grade["false_positive_components"] if seg_grade else None, *b.seg_false_positive_max, "max"),
    }
    g["detection_graded"] = detection_graded
    g["bands"] = bands.__dict__
    with open(cal_dir / "grade.json", "w") as fh:
        json.dump(scrub(g), fh, indent=2, allow_nan=False)
    return g


def scrub(o):
    """JSON-safe copy: numpy scalars to Python, NaN to None, Paths to str."""
    if isinstance(o, dict):
        return {str(k): scrub(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [scrub(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


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

