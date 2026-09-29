"""SEL noise floor: candidates per baseline lesion with NO true expansion, by registration setting and interval.

Uses each subject's calibration synthetic (every lesion unchanged). Registration is redone per setting;
the interval only rescales the Jacobian to %/yr, so it needs no new registration.
Output: study/sel_noise_floor.tsv (one row per subject x setting x interval).
"""
import time
from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from lesiontrack.candidates import sel_candidates
from lesiontrack.config import RegParams, SELParams
from lesiontrack.registration import Timepoint, jacobian_pct_per_year, register_pair
from lesiontrack.tracking import label_lesions

from mscard.spec import read_manifest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "derivatives/mslesseg"
OUT = ROOT / "study"
SETTINGS = {
    "default": RegParams(),
    "smooth": replace(RegParams(), deform_sigma_update="3vox", deform_sigma_total="1vox"),
    "coarse": replace(RegParams(), deform_iterations="50x30x0"),
}
DTS = (0.5, 1.0, 2.0)
man = read_manifest(Path.home() / "repos/lesiontrack/work/manifest_mslesseg.tsv")
rows, tsv = [], OUT / "sel_noise_floor.tsv"
if tsv.exists():
    rows = pd.read_csv(tsv, sep="\t").to_dict("records")
done = {(r["subject"], r["setting"]) for r in rows}
for sub in sorted(man, key=lambda s: int(s[1:])):
    syn = D / sub / "calibration/synthetic"
    bm = D / sub / "measure/brainmask" / f"{man[sub][0].session}_brainmask.nii.gz"
    b0 = man[sub][0]
    if not (syn / "truth.json").exists():
        continue
    brain = b0.brainmask or bm
    bl = Timepoint("BL", b0.t1, b0.flair, b0.mask, 0.0, brain)
    fu = Timepoint("SYN", syn / "T1w.nii.gz", syn / "FLAIR.nii.gz", syn / "mask.nii.gz", 1.0, syn / "brainmask.nii.gz")
    for name, rp in SETTINGS.items():
        if (sub, name) in done:
            continue
        t0 = time.time()
        wd = OUT / "work" / sub / name
        try:
            p = register_pair(bl, fu, wd, replace(rp, threads=4), wd / "greedy.log")
        except Exception as e:  # noqa: BLE001
            print(sub, name, "FAILED", e, flush=True); continue
        m = np.asarray(nib.load(p.baseline_mask).dataobj) > 0
        n_les = int(label_lesions(m, 2).max())
        for dt in DTS:
            exp, _ = jacobian_pct_per_year(p.jacobian, dt)
            c = sel_candidates(exp, m, SELParams())
            n = int(len(np.unique(c)) - 1)
            rows.append({"subject": sub, "setting": name, "dt_years": dt, "baseline_lesions": n_les,
                         "candidates": n, "per_lesion": n / n_les if n_les else None,
                         "lesion_voxels": int(m.sum()), "seconds": round(time.time() - t0)})
        pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
        print(sub, name, rows[-1]["candidates"], f"{time.time()-t0:.0f}s", flush=True)
print("STUDY DONE", flush=True)
