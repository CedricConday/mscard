# mscard

**The longitudinal multiple sclerosis MRI report that grades itself.**

Two (or more) MRI visits in, one page out: new, enlarging, shrinking and resolved
lesions, slowly expanding lesion (SEL) candidates, lesion volume change and brain
volume change. Next to every number, a grade: before the report is written, the
same pipeline is run on the subject's own baseline scan with **known changes
injected** (a known brain contraction, a known set of new lesions, a small
repositioning, gain drift and noise), and the page states how much of that truth
it recovered. A number the pipeline cannot recover on this scan is marked as such.

> **Nothing in this repository is evidence about multiple sclerosis.** It is
> research software. It is not a medical device, it has not been validated
> against clinical outcome, and it makes no diagnostic claim.

Status: first public build, 2026-09-29. Results on the public MSLesSeg cohort
are in the section below as they come in.

## Why

Every clinic that pays for longitudinal MS MRI analysis receives counts of new
and enlarging lesions and a brain volume change, and none of those reports says
how accurate the pipeline was on that scanner, that protocol, that head. The
open pieces exist: [LST-AI](https://github.com/CompImg/LST-AI) segments,
[lesiontrack](https://github.com/CedricConday/lesiontrack) registers, tracks and
finds SEL candidates, [bidsgate](https://github.com/CedricConday/bidsgate)
injects known lesions and atrophy and scores their recovery. mscard is the page
that joins them, and the rule that the calibration is not optional.

## Install

```bash
pip install mscard
```

Python 3.10 or later. Registration comes with lesiontrack (greedy as a wheel,
x86-64 and arm64, CPU only). LST-AI is optional and runs in its container
(`docker` on PATH); without it, give your own lesion masks.

## Use

A manifest, one row per scan, tab separated:

```
subject  session  time_years  t1              flair              mask             brainmask
P1       T1       28.09       P1/T1/T1.nii.gz  P1/T1/FLAIR.nii.gz  P1/T1/MASK.nii.gz  P1/T1/brain.nii.gz
P1       T2       28.51       P1/T2/T1.nii.gz  P1/T2/FLAIR.nii.gz  P1/T2/MASK.nii.gz  P1/T2/brain.nii.gz
```

`time_years` is any monotone clock (age, or years since the first scan).
`mask` and `brainmask` are optional: with `--segmenter lst-ai` the masks come
from LST-AI, and a missing brain mask is estimated (skull-stripped input: the
nonzero support; whole head: bidsgate's morphological estimate). T1w, FLAIR and
masks of one scan must share a grid.

```bash
mscard run manifest.tsv --out derivatives/mscard                 # measure, calibrate, report
mscard run manifest.tsv --out derivatives/mscard --segmenter lst-ai
mscard run manifest.tsv --out derivatives/mscard --no-calibrate  # measure and report only
mscard calibrate manifest.tsv --out derivatives/mscard           # calibration only
mscard report derivatives/mscard                                 # re-render pages
```

Output per subject: `report.html` (self-contained), `measure/mscard.json`
(every number), `measure/lesiontrack/` (masks, label maps, Jacobian, tables),
`calibration/synthetic/` (the synthetic follow-up and its truth),
`calibration/grade.json`. Over the tree: `index.html`, one row per subject.

## What is measured

* **Registration, tracking, SEL candidates**: lesiontrack. Halfway space, greedy
  stationary-velocity deformable registration on T1w and FLAIR jointly, Jacobian
  determinant, lesion matching on the overlap graph; classes follow the
  Vanden Bulcke et al. 2025 bands (enlarging or shrinking beyond ±9 %/yr, stable
  within ±4 %/yr, trend between); SEL candidates by the Elliott et al. 2019 rule.
  A *definite* SEL is a cohort-relative ranking and is not claimed for one subject.
* **Brain volume change**: the mean Jacobian determinant of the same deformation
  field over the baseline brain mask, minus one. When a brain mask exists for
  both visits the mask-volume ratio is shown as a second, ungraded number.

## How the grade is made

For each subject, its baseline scan becomes a synthetic follow-up one modelled
year later (`mscard/calibrate.py`, parameters in `mscard/config.py`):

1. the whole brain contracted by a known factor (default 0.98, i.e. −2 %) with
   bidsgate `inject-atrophy`; every image and mask follows the same field;
2. a known set of new lesions (default 12, 30 to 1500 mm³) with bidsgate
   `inject-lesions`, placed at least 6 mm from every baseline lesion so each is
   truly new;
3. a rigid repositioning (up to 1 degree, 1 voxel), per-image gain (0.9 to 1.1)
   and 2 % noise.

Every baseline lesion is otherwise unchanged. The same pipeline measures
baseline against this synthetic scan and each report line is compared with its
truth:

| Report line | Truth | Grade |
|---|---|---|
| new lesions | the injected set | sensitivity; false new |
| enlarging / shrinking | none | share of unchanged lesions called changed |
| resolved | none | share called resolved |
| SEL candidates | none | candidates per baseline lesion |
| brain volume change | the injected factor | recovery = measured / injected |
| segmentation (LST-AI mode only) | the injected lesions | sensitivity, false positives, Dice |

The traffic-light bands are mscard's own and are printed on every page so the
reader can disagree with them. Green: read the number at face value. Amber: read
it with the calibration table. Red: the pipeline could not recover that kind of
change on this scan; do not use the number.

What the calibration does **not** cover: real new lesions are not ellipsoids and
can sit next to old ones; real atrophy is not uniform; a synthetic year has no
scanner change. It is a floor, not a validation.

## Results

_pending: MSLesSeg cohort, this box, 2026-09-29._

## Citations

* Elliott C, et al. Slowly expanding/evolving lesions as a magnetic resonance imaging marker of chronic active multiple sclerosis lesions. *Mult Scler* 2019.
* Vanden Bulcke C, et al. 2025 (lesion volume change bands, see lesiontrack).
* Wiltgen T, et al. LST-AI: a deep learning ensemble for accurate MS lesion segmentation. *NeuroImage Clin* 2024.
* Guarnera C, et al. MSLesSeg: baseline and benchmarking of a new multiple sclerosis lesion segmentation dataset. 2024/2025.

## Licence

MIT.
