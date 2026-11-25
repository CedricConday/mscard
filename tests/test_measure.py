"""brain_volume_change folds the linear step's volume factor into det J."""

import nibabel as nib
import numpy as np

from mscard.measure import brain_volume_change


def _save(tmp_path, name, arr):
    p = tmp_path / name
    nib.save(nib.Nifti1Image(arr, np.eye(4)), p)
    return p


def test_linear_volume_factor_is_multiplied_in(tmp_path):
    jac = np.full((30, 30, 30), 0.99, np.float32)
    brain = np.zeros((30, 30, 30), np.uint8)
    brain[5:25, 5:25, 5:25] = 1
    j, b = _save(tmp_path, "j.nii.gz", jac), _save(tmp_path, "b.nii.gz", brain)
    only = brain_volume_change(j, b, b, erode_mm=2.0)
    both = brain_volume_change(j, b, b, erode_mm=2.0, linear_volume_factor=0.99)
    assert abs(only["change_pct"] - (-1.0)) < 1e-4
    assert abs(both["change_pct"] - (0.99 * 0.99 - 1) * 100) < 1e-4
    assert abs(both["change_pct_eroded"] - both["change_pct"]) < 1e-4
    assert both["jacobian_only_change_pct"] == only["change_pct"]
