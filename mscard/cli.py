"""mscard command line: run (measure + calibrate + report), calibrate, report."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

from . import __version__
from .config import Bands, Params, SyntheticParams


def _params(a: argparse.Namespace) -> Params:
    sp = SyntheticParams()
    over = {k: getattr(a, k) for k in ("volume_factor", "n_new", "seed", "dt_years") if getattr(a, k, None) is not None}
    if over:
        sp = SyntheticParams(**{**sp.__dict__, **over})
    return Params(synthetic=sp, bands=Bands(), threads=a.threads, engine=a.engine)


def run_subject(subject: str, scans: list, out: Path, segmenter: str, params: Params, do_measure: bool,
                do_calibrate: bool, reference_mask: Path | None = None) -> dict:
    from .calibrate import grade, make_synthetic
    from .measure import measure
    from .segment import ensure_masks, lst_ai
    from .spec import ensure_brainmask

    sub = out / subject
    sub.mkdir(parents=True, exist_ok=True)
    log: list[str] = []
    given_baseline_mask = scans[0].mask
    scans = ensure_masks(scans, segmenter, sub / "segment", params.threads, log)
    status = {"subject": subject}
    if do_measure:
        m = measure(subject, scans, sub / "measure", params.threads, params.engine, log)
        status["measure"] = "ok"
        status["interval_years"] = m["interval_years"]
    if do_calibrate:
        baseline = ensure_brainmask(scans[0], sub / "measure" / "brainmask" / f"{scans[0].session}_brainmask.nii.gz", log)
        cal = sub / "calibration"
        syn, truth = make_synthetic(baseline, cal / "synthetic", params.synthetic, log)
        if segmenter == "lst-ai":
            syn = syn.with_mask(lst_ai(syn, cal / "lst-ai", threads=params.threads))
        m2 = measure(subject, [baseline, syn], cal / "measure", params.threads, params.engine, log)
        ref = given_baseline_mask if (segmenter == "lst-ai" and given_baseline_mask is not None) else None
        g = grade(baseline, syn, cal, m2, truth, segmenter, params.bands, params.threads, reference_mask=ref)
        status["calibrate"] = "ok"
        status["lights"] = g["lights"]
    (sub / "mscard.log").write_text("\n".join(log) + "\n")
    return status


def cmd_run(a: argparse.Namespace, do_measure: bool = True, do_calibrate: bool = True) -> int:
    from .report import gallery, render
    from .spec import read_manifest

    subjects = read_manifest(Path(a.manifest))
    if a.subject:
        want = set(a.subject)
        missing = want - set(subjects)
        if missing:
            raise SystemExit(f"not in manifest: {sorted(missing)}")
        subjects = {k: v for k, v in subjects.items() if k in want}
    if do_calibrate and not do_measure:
        subjects = {k: v[:1] for k, v in subjects.items()}
    out = Path(a.out)
    params = _params(a)
    failed = 0
    for name, scans in subjects.items():
        try:
            st = run_subject(name, scans, out, a.segmenter, params, do_measure, do_calibrate and not a.no_calibrate)
            render(out / name)
            print(f"{name}: {json.dumps(st.get('lights', st))}")
        except Exception as e:  # noqa: BLE001 - one subject must not stop the cohort
            failed += 1
            (out / name).mkdir(parents=True, exist_ok=True)
            (out / name / "FAILED.txt").write_text("".join(traceback.format_exception(e)))
            print(f"{name}: FAILED: {e}", file=sys.stderr)
    gallery(out)
    return 1 if failed else 0


def cmd_report(a: argparse.Namespace) -> int:
    from .report import gallery, render

    out = Path(a.out)
    subs = [p for p in sorted(out.iterdir()) if (p / "measure" / "mscard.json").exists() or (p / "calibration" / "grade.json").exists()]
    if a.subject:
        subs = [p for p in subs if p.name in set(a.subject)]
    for p in subs:
        render(p)
        print(f"{p.name}: {p / 'report.html'}")
    gallery(out)
    print(f"gallery: {out / 'index.html'}")
    return 0


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("manifest", help="TSV: subject, session, time_years, t1, flair [, mask, brainmask, date, scanner]")
    p.add_argument("--out", required=True)
    p.add_argument("--subject", action="append", help="restrict to these subjects (repeatable)")
    p.add_argument("--segmenter", choices=("given", "lst-ai"), default="given",
                   help="'given' uses the manifest's masks; 'lst-ai' runs LST-AI v2 in docker for every scan")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--engine", choices=("greedy", "ants"), default="greedy", help="deformable engine (lesiontrack)")
    p.add_argument("--volume-factor", type=float, dest="volume_factor", help="calibration: injected brain volume factor (default 0.98)")
    p.add_argument("--n-new", type=int, dest="n_new", help="calibration: injected new lesions (default 12)")
    p.add_argument("--seed", type=int, help="calibration seed (default 0), mixed with the subject name")
    p.add_argument("--dt-years", type=float, dest="dt_years", help="calibration: modelled interval (default 1.0)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="mscard", description="The longitudinal MS MRI report that grades itself.")
    ap.add_argument("--version", action="version", version=f"mscard {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="measure every subject, calibrate the pipeline on each baseline, write reports")
    _common(r)
    r.add_argument("--no-calibrate", action="store_true", help="measure and report without the calibration")
    r.set_defaults(func=lambda a: cmd_run(a, True, True))

    c = sub.add_parser("calibrate", help="calibration only: synthetic follow-up of each baseline, grade")
    _common(c)
    c.set_defaults(func=lambda a: cmd_run(a, False, True), no_calibrate=False)

    p = sub.add_parser("report", help="re-render report.html for every subject in an output tree, and the gallery")
    p.add_argument("out")
    p.add_argument("--subject", action="append")
    p.set_defaults(func=cmd_report)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
