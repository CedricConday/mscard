"""Input description: one manifest row per scan, and brain masks when the user has none."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

REQUIRED = ("subject", "session", "time_years", "t1", "flair")
OPTIONAL = ("mask", "brainmask", "date", "scanner")


@dataclass
class Scan:
    subject: str
    session: str
    time_years: float
    t1: Path
    flair: Path
    mask: Path | None = None  # lesion mask; None means mscard must segment
    brainmask: Path | None = None
    date: str = ""
    scanner: str = ""

    def with_mask(self, mask: Path) -> Scan:
        return replace(self, mask=mask)

    def with_brainmask(self, brainmask: Path) -> Scan:
        return replace(self, brainmask=brainmask)


def _p(v: str, base: Path) -> Path | None:
    v = (v or "").strip()
    if not v:
        return None
    p = Path(v)
    return p if p.is_absolute() else (base / p)


def read_manifest(path: Path) -> dict[str, list[Scan]]:
    """TSV with columns subject, session, time_years, t1, flair [, mask, brainmask, date, scanner].

    Relative paths are taken relative to the manifest's directory. Scans are
    sorted by time within a subject; the first is the baseline.
    """
    path = Path(path)
    base = path.parent
    out: dict[str, list[Scan]] = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED if c not in (rd.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: manifest lacks columns {missing}; has {rd.fieldnames}")
        for row in rd:
            if not (row.get("subject") or "").strip():
                continue
            s = Scan(
                subject=row["subject"].strip(),
                session=row["session"].strip(),
                time_years=float(row["time_years"]),
                t1=_p(row["t1"], base),
                flair=_p(row["flair"], base),
                mask=_p(row.get("mask", ""), base),
                brainmask=_p(row.get("brainmask", ""), base),
                date=(row.get("date") or "").strip(),
                scanner=(row.get("scanner") or "").strip(),
            )
            for what in ("t1", "flair"):
                if not getattr(s, what).exists():
                    raise FileNotFoundError(f"{s.subject}/{s.session}: {what} {getattr(s, what)} does not exist")
            for what in ("mask", "brainmask"):
                q = getattr(s, what)
                if q is not None and not q.exists():
                    raise FileNotFoundError(f"{s.subject}/{s.session}: {what} {q} does not exist")
            out.setdefault(s.subject, []).append(s)
    for k, scans in out.items():
        scans.sort(key=lambda x: x.time_years)
        names = [x.session for x in scans]
        if len(set(names)) != len(names):
            raise ValueError(f"{k}: duplicate session names {names}")
    return out


def same_grid(a: Path, b: Path, atol: float = 1e-3) -> bool:
    ia, ib = nib.load(a), nib.load(b)
    return ia.shape == ib.shape and np.allclose(ia.affine, ib.affine, atol=atol)


def check_scan_grid(s: Scan) -> None:
    for what in ("flair", "mask", "brainmask"):
        q = getattr(s, what)
        if q is not None and not same_grid(s.t1, q):
            raise ValueError(f"{s.subject}/{s.session}: {what} is not on the T1w grid; mscard expects one grid per scan")


def ensure_brainmask(s: Scan, out: Path, log: list[str] | None = None) -> Scan:
    """Return the scan with a brain mask, writing one next to ``out`` when none was given.

    Whole-head T1w: bidsgate's morphological estimate (no atlas). Skull-stripped
    T1w (most of the head is exactly zero): the nonzero support, filled.
    """
    if s.brainmask is not None:
        return s
    out = Path(out)
    if out.exists():
        return s.with_brainmask(out)
    from bidsgate.inject_lesions import estimate_brain

    img = nib.load(s.t1)
    t1 = np.asarray(img.dataobj, dtype=np.float32)
    zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
    vox_ml = float(np.prod(zooms)) / 1000.0
    nonzero = t1 > 0
    method = ""
    brain = None
    stripped = nonzero.mean() < 0.45 and 100 <= nonzero.sum() * vox_ml <= 2200
    if stripped:
        lab, n = ndi.label(nonzero)
        if n:
            sizes = np.bincount(lab.ravel())[1:]
            brain = lab == (int(np.argmax(sizes)) + 1)
            brain = ndi.binary_fill_holes(brain)
            method = "nonzero support of a skull-stripped T1w, largest piece, holes filled"
    if brain is None:
        try:
            brain = estimate_brain(t1, zooms, affine=img.affine)
            method = "bidsgate.estimate_brain (morphological, no atlas)"
        except Exception as e:  # noqa: BLE001 - the fallback is stated in the record
            brain = ndi.binary_fill_holes(nonzero)
            method = f"nonzero support (bidsgate.estimate_brain failed: {e})"
    out.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(brain.astype(np.uint8), img.affine), out)
    rec = {"method": method, "volume_ml": float(brain.sum() * vox_ml), "source": str(s.t1)}
    out.with_suffix("").with_suffix(".json").write_text(json.dumps(rec, indent=2))
    if log is not None:
        log.append(f"{s.subject}/{s.session}: brain mask {method}, {rec['volume_ml']:.0f} ml")
    return s.with_brainmask(out)


def to_lesiontrack(s: Scan):
    from lesiontrack.registration import Timepoint

    if s.mask is None:
        raise ValueError(f"{s.subject}/{s.session}: no lesion mask; run the segmenter first")
    return Timepoint(name=s.session, t1=s.t1, flair=s.flair, mask=s.mask, time_years=s.time_years,
                     brainmask=s.brainmask)
