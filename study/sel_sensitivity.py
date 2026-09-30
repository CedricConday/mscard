"""SEL sensitivity to injected expansion, per registration setting (companion to sel_noise_floor.py).

For each MSLesSeg baseline with at least 8 lesions of 30 voxels or more, lesiontrack.synth writes one
synthetic follow-up one year later in which 8 lesions expand (volume x1.15, 1.25, 1.40) and the rest do not.
The pair is registered under the same three settings as the noise-floor study and scored with
lesiontrack.synth.score: sensitivity (expanded lesions with a candidate), false-positive rate (unexpanded
lesions with a candidate) and recovery (measured / injected expansion). Output: study/sel_sensitivity.tsv.
"""

import time
from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from lesiontrack.config import Params, RegParams
from lesiontrack.pipeline import run_subject
from lesiontrack.registration import Timepoint, run_greedy
from lesiontrack.synth import make_followup, score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "study"
SETTINGS = {"default": RegParams(), "smooth": replace(RegParams(), deform_sigma_update="3vox", deform_sigma_total="1vox"),
            "coarse": replace(RegParams(), deform_iterations="50x30x0")}
SRC = Path.home() / "repos/lesiontrack/sourcedata/mslesseg/train"
tsv = OUT / "sel_sensitivity.tsv"
rows = pd.read_csv(tsv, sep="\t").to_dict("records") if tsv.exists() else []
done = {(r["subject"], r["setting"]) for r in rows}
subjects = sorted([p.name for p in SRC.iterdir()], key=lambda s: int(s[1:]))[:30]
for sub in subjects:
    b = SRC / sub / "T1"
    t1, fl, mk, bm = (b / f"{sub}_T1_{x}.nii.gz" for x in ("T1", "FLAIR", "MASK", "brainmask"))
    sdir = OUT / "work_sens" / sub
    try:
        rec = make_followup(t1, fl, mk, sdir / "synthetic", n_expand=8, seed=0, time_fraction=1.0)
    except ValueError as e:
        print(sub, "skipped:", e, flush=True)
        continue
    syn = sdir / "synthetic"
    for name, rp in SETTINGS.items():
        if (sub, name) in done:
            continue
        t0 = time.time()
        tps = [Timepoint("baseline", t1, fl, mk, 0.0, bm),
               Timepoint("synthetic", syn / "synth_T1w.nii.gz", syn / "synth_FLAIR.nii.gz", syn / "synth_mask.nii.gz", 1.0)]
        res = run_subject(sub, tps, sdir / name, Params(reg=replace(rp, threads=4)))
        pair = res.pairs[-1]
        rl = sdir / name / "truth_labels_halfway.nii.gz"
        run_greedy(f"-d 3 -rf {pair.baseline_t1} -ri NN -rt int -rm {syn / 'baseline_lesion_labels.nii.gz'} {rl} -r {pair.halfway_matrix},-1")
        truth = np.asarray(nib.load(rl).dataobj).astype(np.int32)
        exp = np.asarray(nib.load(sdir / name / "expansion_synthetic_pct_per_year.nii.gz").dataobj)
        cand = np.asarray(nib.load(sdir / name / "sel_candidates.nii.gz").dataobj)
        _, m = score(rec, truth, exp, cand, 1.0)
        rows.append({"subject": sub, "setting": name, "sensitivity": m["sensitivity"], "false_positive_rate": m["false_positive_rate"],
                     "recovery_fraction": m["recovery_fraction"], "seconds": round(time.time() - t0)})
        pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
        print(sub, name, round(m["sensitivity"], 2), round(m["false_positive_rate"], 2), round(m["recovery_fraction"], 2), flush=True)
print("STUDY DONE", flush=True)
