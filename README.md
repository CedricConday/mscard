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
(24 subjects) are in the Results section; the gallery of reports is at
https://cedricconday.github.io/mscard/.

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

**MSLesSeg cohort, 2026-09-29, this box (4 CPU threads, greedy engine), expert lesion masks supplied
(`--segmenter given`), calibration factor 0.98, 12 injected lesions per subject.** Every report is
in the gallery: https://cedricconday.github.io/mscard/ (one self-contained page per subject).

| Subject | Years | Baseline lesions | New | Enlarging | Shrinking | Resolved | SEL cand. | Lesion vol %/yr | Brain vol %/yr (eroded) | Boundary | Injected found | False new | Unchanged called changed | SEL/lesion on unchanged | Atrophy recovery |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| P1 | 1.40 | 18 | 4 | 3 | 4 | 7 | 20 | -18.0 | 1.30 (0.01) | amber | 12/12 | 0 | 0/18 | 1.33 | 1.15 |
| P2 | 4.59 | 15 | 2 | 6 | 1 | 0 | 1 | 6.5 | -0.77 (-0.73) | green | 12/12 | 0 | 0/15 | 1.07 | 1.07 |
| P3 | 4.41 | 41 | 5 | 4 | 3 | 21 | 3 | -8.3 | -0.26 (-0.39) | green | 12/12 | 0 | 0/37 | 0.34 | 1.26 |
| P4 | 2.19 | 19 | 4 | 2 | 4 | 5 | 77 | 6.6 | -0.90 (-1.39) | green | 12/12 | 0 | 0/18 | 7.21 | 1.26 |
| P5 | 1.59 | 22 | 4 | 5 | 3 | 9 | 4 | -6.2 | -0.15 (-0.73) | amber | 12/12 | 0 | 0/18 | 0.73 | 1.20 |
| P6 | 2.01 | 33 | 18 | 9 | 6 | 13 | 40 | 14.8 | -2.38 (-2.33) | green | 12/12 | 0 | 0/32 | 1.85 | 1.17 |
| P7 | 0.97 | 21 | 11 | 4 | 3 | 8 | 40 | 30.5 | 0.08 (-0.21) | green | 12/12 | 0 | 0/20 | 3.67 | 1.23 |
| P8 | 1.02 | 17 | 7 | 4 | 7 | 4 | 15 | 10.9 | -3.47 (-2.53) | amber | 12/12 | 0 | 0/16 | 0.41 | 1.15 |
| P9 | 2.15 | 14 | 3 | 4 | 7 | 2 | 4 | 0.6 | -2.46 (-1.23) | amber | 12/12 | 0 | 0/12 | 1.21 | 1.23 |
| P10 | 1.00 | 34 | 11 | 9 | 7 | 12 | 20 | 13.9 | -3.33 (-1.90) | amber | 12/12 | 0 | 0/33 | 0.41 | 1.13 |
| P11 | 0.58 | 14 | 1 | 3 | 4 | 4 | 11 | -44.7 | -3.64 (-5.69) | red | 12/12 | 0 | 0/13 | 0.79 | 1.20 |
| P12 | 3.16 | 32 | 13 | 5 | 8 | 9 | 39 | 5.7 | -1.84 (-1.87) | green | 12/12 | 0 | 0/28 | 2.88 | 1.13 |
| P13 | 1.39 | 38 | 22 | 14 | 14 | 5 | 40 | 4.4 | -2.16 (-2.14) | green | 12/12 | 0 | 0/37 | 1.16 | 1.17 |
| P14 | 3.75 | 16 | 17 | 1 | 2 | 4 | 11 | -0.8 | 0.39 (-0.38) | amber | 12/12 | 0 | 0/16 | 1.31 | 1.10 |
| P19 | 6.03 | 58 | 8 | 3 | 0 | 36 | 19 | 2.5 | 0.54 (0.20) | green | 12/12 | 0 | 1/56 | 1.35 | 1.19 |
| P20 | 4.37 | 8 | 10 | 0 | 1 | 7 | 3 | 24.5 | 0.65 (1.19) | amber | 12/12 | 0 | 0/8 | 0.12 | 1.18 |
| P22 | 0.59 | 21 | 3 | 7 | 9 | 4 | 4 | -14.2 | 0.69 (-1.41) | red | 12/12 | 0 | 0/19 | 0.38 | 1.17 |
| P28 | 0.97 | 8 | 14 | 1 | 3 | 2 | 2 | 25.2 | -8.18 (-6.25) | red | 12/12 | 0 | 0/8 | 0.00 | 1.18 |
| P31 | 1.48 | 28 | 7 | 15 | 4 | 2 | 12 | 34.8 | -0.58 (0.36) | amber | 12/12 | 0 | 0/24 | 0.85 | 1.20 |
| P33 | 2.42 | 41 | 2 | 3 | 8 | 21 | 9 | -16.3 | -0.81 (-0.54) | green | 12/12 | 0 | 0/37 | 1.00 | 1.31 |
| P50 | 1.00 | 47 | 8 | 14 | 10 | 13 | 36 | 298.5 | -0.37 (0.42) | amber | 12/12 | 0 | 0/43 | 1.11 | 1.46 |
| P51 | 3.01 | 43 | 18 | 6 | 2 | 26 | 33 | 9.5 | -1.29 (-1.72) | green | 12/12 | 0 | 0/40 | 1.70 | 1.34 |
| P52 | 0.16 | 26 | 3 | 5 | 4 | 15 | 24 | -87.5 | 29.48 (16.12) | red | 12/12 | 0 | 0/25 | 1.00 | 1.25 |
| P53 | 1.17 | 22 | 21 | 0 | 1 | 21 | 21 | -7.6 | -0.97 (-1.41) | green | 12/12 | 0 | 0/21 | 1.14 | 1.13 |

Subjects: 24. Injected new lesions reported: 288 of 288; false new: 0. Unchanged lesions called enlarging or shrinking: 1 of 594; called resolved: 0. SEL candidates per unchanged lesion: median 1.09 (range 0.00 to 7.21). Atrophy recovery: median 1.18 (range 1.07 to 1.46). Boundary check on the real pairs: {'green': 11, 'amber': 9, 'red': 4}.

P49 is not in the table: its released follow-up T1w is empty (lesiontrack refuses it; reported upstream
as MSLesSeg-2024 issue #1). P50's baseline leaves too little white matter for lesion placement with a 6 mm
margin around its existing lesions; since 2026-09-30 the calibration steps the margin down and records it
(P50: 4 mm).

What the numbers say, and what they do not:

* **New lesions, supplied masks: matching only.** With masks drawn by an annotator, the synthetic
  follow-up mask contains the injected lesions by construction, so 276 of 276 says that
  lesiontrack's tracking calls an injected lesion "new" rather than a fragment or a merge, and
  nothing about detection. That is why the badge is grey. The LST-AI mode grades detection; its
  numbers follow below as they land.
* **Enlarging, shrinking, resolved: 1 false call in 594 lesions, 0 false resolved.** Each baseline
  lesion is judged against its own truth volume in the synthetic mask, with lesiontrack's ±9 %/yr
  and ±4 %/yr bands. Registration and tracking do not invent change on a one-year synthetic
  interval. The real pairs still show many resolved and shrinking lesions (P1: 7 resolved of 18);
  the calibration cannot tell annotation drift between visits from biology, and the page says so.
* **Brain volume change: recovery 1.07 to 1.46, median 1.18.** The deformation field over-reads a
  uniform 2 % contraction by about a fifth, consistently (lesiontrack's lesion-expansion backtest
  gave 1.16). On the real pairs the boundary check is green for 11, amber for 9, red for 4: where it
  is red, the full-mask and eroded-mask numbers disagree by more than 1.5 points per year and the
  page calls the number rim-dominated.
* **SEL candidates: red for 18 of 24.** With no lesion expanding, the Elliott et al. 2019 candidate
  rule still marks a median 1.07 candidates per baseline lesion. That is the noise floor of the
  candidate stage at this registration; the definite/possible split is a cohort ranking on top of
  it. This is the finding a reader of any SEL count should have next to it.

### LST-AI mode (mscard runs the segmenter): two subjects, 2026-09-29

Same subjects, same synthetic follow-ups, but the lesion masks of every scan (real and synthetic)
come from LST-AI v2.0.0rc1 (CPU, fast mode), so detection is graded too. Reports: P2 and P20 in the
gallery under `lstai/`.

| Subject | Injected found by LST-AI | False-positive components | Dice on injected | Baseline Dice vs expert mask | Reported new | False new | Baseline lesions falsely called changed | Falsely resolved | Atrophy recovery |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| P2 | 8/12 | 7 | 0.52 | 0.73 | 8/12 | 3 | 6/15 | 0 | 1.07 |
| P20 | 2/12 | 2 | 0.03 | 0.15 | 2/12 | 2 | 0/1 | 0 | 1.17 |

This is the page grading itself. With supplied masks the tracking lines were green; with LST-AI
the new-lesion badge is red on both subjects, the false-change badge is red on P2 (segmentation
differences between two visits read as lesion change), and on P20 LST-AI agrees with the expert
mask at Dice 0.15 at baseline, so nothing on that page should be read at face value, and the page
says so. The synthetic lesions are ellipsoids with a fixed FLAIR contrast, so this is a floor for
the segmenter, not its accuracy on real lesions. A third subject (P3) was stopped before its
synthetic scan.

## FHIR export

`mscard fhir derivatives/mscard/P1 --patient Patient/123 --out P1_bundle.json` writes the report as a FHIR R4B Bundle:
one DiagnosticReport and one Observation per finding, each carrying the calibration light as its `interpretation`,
so a hospital system receives the number together with whether the same pipeline recovered that kind of change on
this patient's own scan. Validated against the R4B models of `fhir.resources` (`pip install mscard[fhir]`); status
`preliminary`; finding codes are local because no LOINC code exists for them. (Formerly the separate package
mscard-fhir.)

## Citations

* Elliott C, et al. Slowly expanding/evolving lesions as a magnetic resonance imaging marker of chronic active multiple sclerosis lesions. *Mult Scler* 2019.
* Vanden Bulcke C, et al. 2025 (lesion volume change bands, see lesiontrack).
* Wiltgen T, et al. LST-AI: a deep learning ensemble for accurate MS lesion segmentation. *NeuroImage Clin* 2024.
* Guarnera C, et al. MSLesSeg: baseline and benchmarking of a new multiple sclerosis lesion segmentation dataset. 2024/2025.

## Licence

MIT.
