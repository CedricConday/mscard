"""Phantom tests: no data download, a synthetic head, every stage exercised except the segmenter."""

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from mscard.calibrate import light
from mscard.config import Bands, SyntheticParams
from mscard.spec import Scan, ensure_brainmask, read_manifest


def _phantom(tmp: Path, shape=(96, 112, 96)) -> Scan:
    rng = np.random.default_rng(0)
    zz, yy, xx = np.meshgrid(*[np.arange(n) for n in shape], indexing="ij")
    c = np.array(shape) / 2
    r = np.sqrt(((zz - c[0]) / 40) ** 2 + ((yy - c[1]) / 48) ** 2 + ((xx - c[2]) / 40) ** 2)
    brain = r < 1
    t1 = np.where(brain, 200.0, 0.0).astype(np.float32)
    t1[brain] += rng.normal(0, 5, int(brain.sum()))
    fl = np.where(brain, 80.0, 0.0).astype(np.float32)
    mask = np.zeros(shape, np.uint8)
    mask[40:46, 52:58, 40:46] = 1
    t1[mask > 0] -= 40
    fl[mask > 0] += 60
    aff = np.diag([1.0, 1.0, 1.0, 1.0])
    for name, arr in (("T1w", t1), ("FLAIR", fl), ("mask", mask)):
        nib.save(nib.Nifti1Image(arr, aff), tmp / f"{name}.nii.gz")
    return Scan("phantom", "bl", 0.0, tmp / "T1w.nii.gz", tmp / "FLAIR.nii.gz", tmp / "mask.nii.gz")


def test_manifest_roundtrip(tmp_path):
    s = _phantom(tmp_path)
    (tmp_path / "m.tsv").write_text(
        "subject\tsession\ttime_years\tt1\tflair\tmask\tbrainmask\n"
        f"phantom\tfu\t1.0\t{s.t1}\t{s.flair}\t{s.mask}\t\n"
        f"phantom\tbl\t0.0\t{s.t1.name}\t{s.flair.name}\t{s.mask.name}\t\n"
    )
    m = read_manifest(tmp_path / "m.tsv")
    assert list(m) == ["phantom"]
    assert [x.session for x in m["phantom"]] == ["bl", "fu"]  # sorted by time
    assert m["phantom"][0].t1 == s.t1  # relative path resolved against the manifest
    assert m["phantom"][0].brainmask is None


def test_manifest_rejects_missing_file(tmp_path):
    (tmp_path / "m.tsv").write_text("subject\tsession\ttime_years\tt1\tflair\nA\tbl\t0\tnope.nii.gz\tnope.nii.gz\n")
    with pytest.raises(FileNotFoundError):
        read_manifest(tmp_path / "m.tsv")


def test_brainmask_from_skull_stripped(tmp_path):
    s = _phantom(tmp_path)
    log = []
    s2 = ensure_brainmask(s, tmp_path / "bm.nii.gz", log)
    bm = np.asarray(nib.load(s2.brainmask).dataobj) > 0
    t1 = np.asarray(nib.load(s.t1).dataobj) > 0
    assert bm.sum() == t1.sum() and log and "skull-stripped" in log[0]


def test_synthetic_truth_is_consistent(tmp_path):
    from mscard.calibrate import make_synthetic

    s = ensure_brainmask(_phantom(tmp_path), tmp_path / "bm.nii.gz")
    sp = SyntheticParams(n_new=3, volume_factor=0.95, rigid_rot_deg=0.5, rigid_trans_vox=0.5, exclusion_mm=3.0)
    syn, truth = make_synthetic(s, tmp_path / "syn", sp)
    assert truth["n_new"] == 3 and len(truth["new_lesions"]) == 3
    new = np.asarray(nib.load(tmp_path / "syn" / "newlesions_truth.nii.gz").dataobj)
    mask = np.asarray(nib.load(syn.mask).dataobj) > 0
    assert set(np.unique(new)) == {0, 1, 2, 3}
    assert mask[new > 0].all()  # the truth follow-up mask contains every injected lesion
    # the contracted brain is smaller by about the factor
    b0 = (np.asarray(nib.load(s.brainmask).dataobj) > 0).sum()
    b1 = (np.asarray(nib.load(syn.brainmask).dataobj) > 0).sum()
    assert abs(b1 / b0 - 0.95) < 0.03
    # the cache returns the same thing without rewriting
    syn2, truth2 = make_synthetic(s, tmp_path / "syn", sp)
    assert truth2["seed"] == truth["seed"]


def test_lights():
    b = Bands()
    assert light(0.95, *b.new_sensitivity_min, "min") == "green"
    assert light(0.75, *b.new_sensitivity_min, "min") == "amber"
    assert light(0.2, *b.new_sensitivity_min, "min") == "red"
    assert light(0, *b.false_new_max, "max") == "green"
    assert light(1.1, *b.atrophy_recovery_range, "range") == "green"
    assert light(1.4, *b.atrophy_recovery_range, "range") == "amber"
    assert light(None, *b.atrophy_recovery_range, "range") == "none"
    assert light(float("nan"), *b.false_new_max, "max") == "none"
