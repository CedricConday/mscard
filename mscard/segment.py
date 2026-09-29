"""Lesion segmentation. Either the user supplies masks, or mscard runs LST-AI in its container."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

import nibabel as nib
import numpy as np

from .config import LST_AI_IMAGE
from .spec import Scan, same_grid

LST_OUT_NAME = "space-flair_seg-lst.nii.gz"


def docker_available() -> bool:
    return shutil.which("docker") is not None


def lst_ai(s: Scan, out_dir: Path, threads: int = 4, fast: bool = True, image: str = LST_AI_IMAGE) -> Path:
    """Run LST-AI v2 (segment only) on one scan; return the binary lesion mask on the scan's grid.

    The container runs as root and keeps HD-BET weights in /root, so the output
    directory is chowned back afterwards through the same image. The result is
    cached: an existing mask is returned without running anything.
    """
    out_dir = Path(out_dir)
    mask = out_dir / f"{s.subject}_{s.session}_desc-lstai_mask.nii.gz"
    if mask.exists():
        return mask
    if not docker_available():
        raise RuntimeError("docker is not on PATH; give lesion masks in the manifest or install docker")
    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    name = f"mscard-lst-{s.subject}-{s.session}-{os.getpid()}"
    subprocess.run(["docker", "rm", "-f", name], capture_output=True, check=False)
    t0 = time.time()
    cmd = [
        "docker", "run", "--name", name,
        "-v", f"{s.t1.parent.resolve()}:/in1:ro",
        "-v", f"{s.flair.parent.resolve()}:/in2:ro",
        "-v", f"{out_dir.resolve()}:/out",
        "-v", f"{work.resolve()}:/tmpd",
        image,
        "--t1", f"/in1/{s.t1.name}", "--flair", f"/in2/{s.flair.name}",
        "--output", "/out", "--temp", "/tmpd", "--device", "cpu", "--segment_only", "--threads", str(threads),
    ]
    if fast:
        cmd.append("--fast-mode")
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (out_dir / "lst-ai.log").write_text(res.stdout + "\n" + res.stderr)
    subprocess.run(["docker", "rm", name], capture_output=True, check=False)
    subprocess.run(["docker", "run", "--rm", "-v", f"{out_dir.resolve()}:/p", "--entrypoint", "chown", image,
                    "-R", f"{os.getuid()}:{os.getgid()}", "/p"], capture_output=True, check=False)
    seg = out_dir / LST_OUT_NAME
    if res.returncode != 0 or not seg.exists():
        raise RuntimeError(f"LST-AI failed for {s.subject}/{s.session} (exit {res.returncode}); see {out_dir / 'lst-ai.log'}")
    if not same_grid(s.t1, seg):
        raise RuntimeError(f"LST-AI wrote {seg} on a grid other than the T1w's; mscard expects FLAIR and T1w on one grid")
    img = nib.load(seg)
    m = (np.asarray(img.dataobj) > 0).astype(np.uint8)
    nib.save(nib.Nifti1Image(m, img.affine), mask)
    (out_dir / "lst-ai.json").write_text(
        f'{{"image": "{image}", "fast_mode": {str(fast).lower()}, "threads": {threads}, '
        f'"seconds": {time.time() - t0:.0f}, "lesion_voxels": {int(m.sum())}}}\n'
    )
    shutil.rmtree(work, ignore_errors=True)
    return mask


def ensure_masks(scans: list[Scan], segmenter: str, out_dir: Path, threads: int, log: list[str]) -> list[Scan]:
    """Give every scan a lesion mask according to the segmenter choice."""
    done = []
    for s in scans:
        if segmenter == "given":
            if s.mask is None:
                raise ValueError(f"{s.subject}/{s.session}: segmenter 'given' but the manifest has no mask")
            done.append(s)
        elif segmenter == "lst-ai":
            m = lst_ai(s, Path(out_dir) / "lst-ai" / s.session, threads=threads)
            log.append(f"{s.subject}/{s.session}: LST-AI mask {m}")
            done.append(s.with_mask(m))
        else:
            raise ValueError(f"unknown segmenter {segmenter!r}; use 'given' or 'lst-ai'")
    return done
