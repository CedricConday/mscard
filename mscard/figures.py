"""Overview figure per pair: baseline, warped follow-up, expansion map. Each image windowed on its own."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load(p: Path) -> np.ndarray:
    return np.asarray(nib.load(p).dataobj)


def _window(img: np.ndarray) -> tuple[float, float]:
    nz = img[img > 0]
    return (0.0, float(np.percentile(nz, 99))) if nz.size else (0.0, 1.0)


def overview(lt_dir: Path, follow_up: str, out: Path, n_slices: int = 4) -> Path:
    """Axial slices through the most lesioned levels, at least a quarter of the span apart."""
    lt_dir = Path(lt_dir)
    t1 = _load(lt_dir / "reg" / f"{follow_up}_baseline_T1w_halfway.nii.gz").astype(np.float32)
    fu = _load(lt_dir / "reg" / f"{follow_up}_followup_T1w_warped.nii.gz").astype(np.float32)
    exp = _load(lt_dir / f"expansion_{follow_up}_pct_per_year.nii.gz").astype(np.float32)
    les = _load(lt_dir / "baseline_lesion_labels.nii.gz")
    cand = _load(lt_dir / "sel_candidates.nii.gz")
    new = _load(lt_dir / f"new_or_enlarging_{follow_up}.nii.gz")
    per_slice = (les > 0).sum(axis=(0, 1))
    min_gap = max(3, les.shape[2] // (4 * n_slices))
    slices: list[int] = []
    for z in np.argsort(per_slice)[::-1]:
        if per_slice[z] == 0:
            break
        if all(abs(int(z) - s) >= min_gap for s in slices):
            slices.append(int(z))
        if len(slices) == n_slices:
            break
    if not slices:
        slices = [les.shape[2] // 2]
    slices.sort()
    w_t1, w_fu = _window(t1), _window(fu)
    fig, axes = plt.subplots(len(slices), 3, figsize=(10, 3.3 * len(slices)))
    axes = np.atleast_2d(axes)
    for row, z in zip(axes, slices):
        for ax, img, w, title in zip(row, (t1, fu, t1), (w_t1, w_fu, w_t1),
                                     ("baseline T1w", f"{follow_up} T1w warped to baseline", "expansion in lesions (%/yr)")):
            ax.imshow(np.rot90(img[:, :, z]), cmap="gray", vmin=w[0], vmax=w[1])
            ax.set_title(f"{title}  z={z}", fontsize=8)
            ax.axis("off")
        row[0].contour(np.rot90(les[:, :, z] > 0), levels=[0.5], colors="cyan", linewidths=0.6)
        row[1].contour(np.rot90(new[:, :, z] > 0), levels=[0.5], colors="yellow", linewidths=0.6)
        heat = np.rot90(np.where(les[:, :, z] > 0, exp[:, :, z], np.nan))
        row[2].imshow(heat, cmap="coolwarm", vmin=-25, vmax=25, alpha=0.9)
        row[2].contour(np.rot90(cand[:, :, z] > 0), levels=[0.5], colors="lime", linewidths=0.8)
    fig.suptitle("cyan: baseline lesions   yellow: new or enlarging voxels   lime: SEL candidates", fontsize=9)
    fig.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out
