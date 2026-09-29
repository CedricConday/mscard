"""One self-contained HTML report per subject, and a gallery over an output tree."""

from __future__ import annotations

import base64
import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

COLOURS = {"green": "#2e7d32", "amber": "#ef6c00", "red": "#c62828", "none": "#757575"}
WORDS = {"green": "read at face value", "amber": "read with the calibration", "red": "not trustworthy here", "none": "not graded"}

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1080px;margin:24px auto;padding:0 16px;color:#1a1a1a;line-height:1.45}
h1{font-size:1.6rem;margin:0 0 4px} h2{font-size:1.15rem;margin:28px 0 8px;border-bottom:1px solid #ddd;padding-bottom:4px}
.meta{color:#555;font-size:.92rem} .cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:16px 0}
.card{border:1px solid #ddd;border-radius:8px;padding:10px 12px;background:#fafafa} .card .n{font-size:1.7rem;font-weight:600;line-height:1.1}
.card .l{font-size:.85rem;color:#444} .badge{display:inline-block;border-radius:4px;color:#fff;font-size:.72rem;padding:1px 6px;margin-top:6px;letter-spacing:.02em}
table{border-collapse:collapse;width:100%;font-size:.9rem} th,td{border-bottom:1px solid #e5e5e5;padding:5px 8px;text-align:left;vertical-align:top}
th{background:#f3f3f3} td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.box{border-left:4px solid #999;background:#f7f7f7;padding:10px 14px;margin:12px 0;font-size:.92rem}
.warn{border-left-color:#c62828;background:#fdf3f3} img{max-width:100%;border:1px solid #ddd;border-radius:6px}
.small{font-size:.82rem;color:#555} .dot{display:inline-block;width:.8em;height:.8em;border-radius:50%;vertical-align:-1px;margin-right:6px}
"""


def _f(x, d=1, suffix="") -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    return f"{x:.{d}f}{suffix}"


def _i(x) -> str:
    return "–" if x is None else f"{int(x)}"


def _badge(light: str, text: str | None = None) -> str:
    return f'<span class="badge" style="background:{COLOURS[light]}">{html.escape(text or WORDS[light])}</span>'


def _dot(light: str) -> str:
    return f'<span class="dot" style="background:{COLOURS.get(light, COLOURS["none"])}"></span>'


def _img(path: Path | None) -> str:
    if not path or not Path(path).exists():
        return '<p class="small">image not available</p>'
    b = base64.b64encode(Path(path).read_bytes()).decode()
    return f'<img src="data:image/png;base64,{b}" alt="{html.escape(Path(path).name)}">'


def _load(p: Path) -> dict | None:
    return json.loads(p.read_text()) if p.exists() else None


def render(subject_dir: Path) -> Path:
    subject_dir = Path(subject_dir)
    m = _load(subject_dir / "measure" / "mscard.json")
    g = _load(subject_dir / "calibration" / "grade.json")
    t = _load(subject_dir / "calibration" / "synthetic" / "truth.json")
    m2 = _load(subject_dir / "calibration" / "measure" / "mscard.json")
    failed = subject_dir / "FAILED.txt"
    subject = subject_dir.name
    L = (g or {}).get("lights", {})
    parts = [(f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
              f"<title>mscard {html.escape(subject)}</title><style>{CSS}</style></head><body>")]
    parts.append(f"<h1>mscard report: {html.escape(subject)}</h1>")
    if m:
        s0, s1 = m["scans"][0], m["scans"][-1]
        parts.append(f"<div class='meta'>Baseline {html.escape(s0['session'])}{(' (' + html.escape(s0['date']) + ')') if s0['date'] else ''}"
                     f" to follow-up {html.escape(s1['session'])}{(' (' + html.escape(s1['date']) + ')') if s1['date'] else ''}, "
                     f"interval {_f(m['interval_years'], 2)} years"
                     f"{', ' + html.escape(s0['scanner']) if s0.get('scanner') else ''}. "
                     f"Lesion masks: {'supplied' if (g or {}).get('segmenter', 'given') == 'given' else 'LST-AI v2, run by mscard'}. "
                     f"mscard {html.escape(m['mscard_version'])}, lesiontrack {html.escape(m['lesiontrack_version'])}, "
                     f"generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}.</div>")
    if failed.exists():
        parts.append(f"<div class='box warn'><b>This subject failed.</b><pre>{html.escape(failed.read_text()[-2000:])}</pre></div>")
    parts.append("<div class='box'><b>Research software, not a medical device.</b> Nothing on this page is a diagnosis. "
                 "Every number below is produced by an open pipeline and, where a coloured badge is shown, the same pipeline was "
                 "run on this subject's own baseline scan with known changes injected, and the badge says how well it recovered them. "
                 "A green badge means the number can be read at face value on this scan; amber, read it together with the calibration "
                 "table; red, the pipeline could not recover that kind of change here and the number should not be used.</div>")

    if m:
        les, br = m["lesions"], m["brain"]
        parts.append("<h2>Findings</h2><div class='cards'>")
        parts.append(_card(_i(les["new"]), "new lesions" + (f" (+{les['adjacent_fragment']} adjacent fragments)" if les["adjacent_fragment"] else ""),
                           L.get("new_sensitivity", "none"), _light_text(g, "new")))
        parts.append(_card(_i(les["enlarging"]), "enlarging lesions", L.get("false_change_rate", "none"), _light_text(g, "change")))
        parts.append(_card(_i(les["shrinking"]), "shrinking lesions", L.get("false_change_rate", "none"), _light_text(g, "change")))
        parts.append(_card(_i(les["resolved"]), "resolved lesions", L.get("false_resolved_rate", "none"), _light_text(g, "resolved")))
        parts.append(_card(_f(les["volume_change_pct_per_year"], 1, " %/yr"), f"lesion volume change ({_f(les['baseline_volume_mm3'] / 1000, 2)} to {_f(les['followup_volume_mm3'] / 1000, 2)} ml)",
                           L.get("false_change_rate", "none"), _light_text(g, "change")))
        parts.append(_card(_i(m["sel"]["candidates"]), "slowly expanding lesion candidates", L.get("sel_candidates_per_lesion", "none"), _light_text(g, "sel")))
        bl = br["jacobian"].get("boundary_light", "none")
        al = L.get("atrophy_recovery", "none")
        worst = max((al, bl), key=lambda x: ["none", "green", "amber", "red"].index(x))
        parts.append(_card(_f(br["jacobian"]["change_pct_per_year"], 2, " %/yr"), "brain volume change (deformation field)",
                           worst, _light_text(g, "atrophy") + f"; boundary check: {_f(br['jacobian'].get('change_pct_per_year_eroded'), 2, ' %/yr')} with the mask eroded {_f(br['jacobian'].get('erode_mm'), 0, ' mm')}"))
        if br.get("mask_ratio"):
            parts.append(_card(_f(br["mask_ratio"]["change_pct_per_year"], 2, " %/yr"), "brain volume change (mask ratio)", "none", "not graded"))
        parts.append("</div>")
        parts.append(f"<p class='small'>Brain volume change is the mean Jacobian determinant over the baseline brain mask. Its boundary check repeats the mean with the mask eroded "
                     f"{_f(br['jacobian'].get('erode_mm'), 0)} mm: {_f(br['jacobian'].get('change_pct'), 2, ' %')} full against {_f(br['jacobian'].get('change_pct_eroded'), 2, ' %')} eroded over the interval "
                     f"(disagreement {_f(br['jacobian'].get('boundary_disagreement_pct'), 2, ' points')}: green ≤ 0.5, amber ≤ 1.5, red beyond). A change that lives only in the outer rim is a mask "
                     f"difference between the two visits, which the calibration cannot see because a synthetic follow-up is stripped identically.</p>")
        parts.append(f"<p class='small'>Baseline: {_i(les['baseline_count'])} lesions, {_f(les['baseline_volume_mm3'] / 1000, 2)} ml. "
                     f"Classes follow Vanden Bulcke et al. 2025 bands: enlarging or shrinking beyond ±9 %/yr, stable within ±4 %/yr, "
                     f"trend between ({_i(les['trend_up'])} trend up, {_i(les['trend_down'])} trend down, {_i(les['stable'])} stable). "
                     f"SEL candidates are the Elliott et al. 2019 Jacobian rule; a definite SEL is a cohort-relative ranking and needs a cohort.</p>")

    if g and t:
        parts.append("<h2>How this report was graded</h2>")
        parts.append(f"<p>The baseline scan of this subject was turned into a synthetic follow-up one modelled year later: the whole brain "
                     f"contracted by <b>{_f(t['injected_brain_change_pct'], 1, ' %')}</b>, <b>{t['n_new']} new lesions</b> "
                     f"({_f(min(x['volume_mm3'] for x in t['new_lesions']), 0)} to {_f(max(x['volume_mm3'] for x in t['new_lesions']), 0)} mm³) "
                     f"placed at least {t['params']['exclusion_mm']:.0f} mm from every existing lesion, a rigid repositioning "
                     f"({', '.join(_f(v, 2) for v in t['rigid']['rot_deg'])} degrees; {', '.join(_f(v, 2) for v in t['rigid']['trans_vox'])} voxels), "
                     f"gains {', '.join(_f(v, 2) for v in t['gains'])} and {t['params']['noise_frac'] * 100:.0f} % noise. "
                     f"Every baseline lesion is otherwise unchanged. The same pipeline then measured baseline against this synthetic scan.</p>")
        n, u, b, s = g["new_lesions"], g["unchanged_lesions"], g["brain"], g["sel"]
        rows = [
            ("New lesions", f"{n['injected']} injected", (f"{n['reported']} reported ({n['reported_as_new']} new, {n['reported_as_fragment']} adjacent fragment); "
             f"{n['false_new']} false new, {n['false_fragment']} false fragments; largest missed {_f(n['largest_missed_mm3'], 0, ' mm³')}"),
             L["new_sensitivity"], f"sensitivity {_f(n['sensitivity'], 2)}"),
            ("False new lesions", "0", f"{n['false_new']}", L["false_new"], ""),
            ("Enlarging or shrinking calls on unchanged lesions", f"0 of {u['n']}", (f"{u['called_enlarging_or_shrinking']} ({_f((u['false_change_rate'] or 0) * 100, 1, ' %')}); "
             f"{u['called_trend']} trend calls"), L["false_change_rate"], ""),
            ("Resolved calls on unchanged lesions", f"0 of {u['n']}", f"{u['called_resolved']}", L["false_resolved_rate"], ""),
            ("SEL candidates on unchanged lesions", "0", f"{s['candidates']} ({_f(s['per_baseline_lesion'], 2)} per baseline lesion)", L["sel_candidates_per_lesion"], ""),
            ("Brain volume change", _f(b["injected_change_pct"], 2, " %"), f"{_f(b['jacobian_change_pct'], 2, ' %')} from the deformation field (recovery {_f(b['jacobian_recovery'], 2)})"
             + (f"; {_f(b['mask_ratio_change_pct'], 2, ' %')} from the mask ratio (recovery {_f(b['mask_ratio_recovery'], 2)})" if b.get("mask_ratio_change_pct") is not None else ""),
             L["atrophy_recovery"], ""),
        ]
        if g.get("segmentation"):
            sg = g["segmentation"]
            rows.append(("Segmentation of the injected lesions (LST-AI)", f"{sg['lesions']} lesions", f"{sg['detected']} found, Dice {_f(sg['dice'], 2)}, "
                         f"{sg['false_positive_components']} false-positive components ({_f(sg['false_positive_volume_mm3'], 0, ' mm³')})"
                         + (f"; baseline Dice against the supplied mask {_f(sg.get('baseline_reference_dice'), 2)}" if sg.get("baseline_reference_dice") is not None else ""),
                         L["seg_sensitivity"], f"sensitivity {_f(sg['sensitivity'], 2)}"))
        parts.append("<table><tr><th>Report line</th><th>Truth</th><th>Pipeline reported</th><th>Reading</th></tr>")
        for name, truth, got, light, extra in rows:
            parts.append(f"<tr><td>{html.escape(name)}</td><td>{html.escape(truth)}</td><td>{html.escape(got)}</td>"
                         f"<td>{_dot(light)}{html.escape(WORDS[light])}{(' (' + html.escape(extra) + ')') if extra else ''}</td></tr>")
        parts.append("</table>")
        if n.get("by_size"):
            parts.append("<p class='small'>Injected lesions reported by size: " + "; ".join(
                f"{k} mm³: {v['reported']} of {v['n']}" for k, v in n["by_size"].items()) + ".</p>")
        bands = g.get("bands", {})
        parts.append(f"<p class='small'>Bands (mscard's own, stated so you can disagree): new-lesion sensitivity green ≥ {bands.get('new_sensitivity_min', ['?'])[0]}, "
                     f"false new ≤ {bands.get('false_new_max', ['?'])[0]}, brain-volume recovery {bands.get('atrophy_recovery_range', [['?', '?']])[0]}, "
                     f"false change rate ≤ {bands.get('false_change_rate_max', ['?'])[0]}, false resolved ≤ {bands.get('false_resolved_rate_max', ['?'])[0]}, "
                     f"SEL candidates per lesion ≤ {bands.get('sel_candidates_per_lesion_max', ['?'])[0]}; amber up to the second value of each.</p>")
    elif m:
        parts.append("<h2>How this report was graded</h2><p>No calibration was run for this subject (<code>--no-calibrate</code>), so no badge is coloured.</p>")

    if m:
        parts.append("<h2>Images</h2>")
        parts.append(f"<p class='small'>Real pair, baseline {html.escape(m['baseline'])} to {html.escape(m['follow_up'])}: baseline lesions (cyan), "
                     "new voxels on the warped follow-up (yellow), Jacobian expansion inside lesions with SEL candidate outlines (lime).</p>")
        parts.append(_img(Path(m["images"]["overview"])))
        if m2:
            parts.append("<p class='small'>Calibration pair, baseline to synthetic follow-up, same rendering.</p>")
            parts.append(_img(Path(m2["images"]["overview"])))
        parts.append("<h2>Lesions</h2>")
        rows = [r for r in m["lesions"]["table"] if r["class"] != "stable"][:40]
        parts.append("<table><tr><th>Group</th><th>Class</th><th class='num'>Baseline mm³</th><th class='num'>Follow-up mm³</th><th class='num'>%/yr</th><th>Voxel centre (halfway grid)</th></tr>")
        for r in rows:
            parts.append(f"<tr><td>{r['group_id']}</td><td>{html.escape(str(r['class']))}</td><td class='num'>{_f(r['volume_baseline_mm3'], 0)}</td>"
                         f"<td class='num'>{_f(r['volume_followup_mm3'], 0)}</td><td class='num'>{_f(r['pct_per_year'], 1)}</td>"
                         f"<td>{_f(r['cx'], 0)}, {_f(r['cy'], 0)}, {_f(r['cz'], 0)}</td></tr>")
        parts.append("</table>")
        if len(m["lesions"]["table"]) > len(rows):
            parts.append(f"<p class='small'>{len(m['lesions']['table']) - len(rows)} stable lesions not listed; the full table is lesion_tracking.tsv.</p>")

    parts.append("<h2>Method</h2><ul class='small'>"
                 "<li>Registration, Jacobian, lesion tracking and SEL candidates: lesiontrack (halfway space, greedy stationary-velocity deformable registration on T1w and FLAIR jointly; "
                 "Elliott et al. 2019 SEL rule; Vanden Bulcke et al. 2025 change bands).</li>"
                 "<li>Brain volume change: mean Jacobian determinant of the same deformation field over the baseline brain mask, minus one.</li>"
                 "<li>Calibration: bidsgate inject-atrophy and inject-lesions on the subject's own baseline; lesiontrack's synthetic rigid perturbation; grade in mscard/calibrate.py.</li>"
                 "<li>Segmentation, when mscard ran it: LST-AI v2.0.0rc1 (CPU container, segment only, fast mode).</li>"
                 f"<li>mscard {html.escape(__version__)}: https://github.com/CedricConday/mscard</li></ul>")
    parts.append("</body></html>")
    out = subject_dir / "report.html"
    out.write_text("\n".join(parts))
    return out


def _card(n: str, label: str, light: str, text: str) -> str:
    return f"<div class='card'><div class='n'>{n}</div><div class='l'>{html.escape(label)}</div>{_badge(light, text)}</div>"


def _light_text(g: dict | None, what: str) -> str:
    if not g:
        return "not graded"
    n, u, b, s = g["new_lesions"], g["unchanged_lesions"], g["brain"], g["sel"]
    if what == "new":
        return f"calibration: {n['reported']} of {n['injected']} injected found, {n['false_new']} false"
    if what == "change":
        return f"calibration: {u['called_enlarging_or_shrinking']} of {u['n']} unchanged lesions called changed"
    if what == "resolved":
        return f"calibration: {u['called_resolved']} of {u['n']} unchanged lesions called resolved"
    if what == "sel":
        return f"calibration: {s['candidates']} candidates with no true expansion"
    if what == "atrophy":
        return f"calibration: {_f(b['injected_change_pct'], 1, ' %')} injected, {_f(b['jacobian_change_pct'], 1, ' %')} measured"
    return "not graded"


def gallery(out: Path) -> Path:
    out = Path(out)
    rows = []
    for p in sorted(out.iterdir()):
        if not p.is_dir():
            continue
        m = _load(p / "measure" / "mscard.json")
        g = _load(p / "calibration" / "grade.json")
        if not m and not g and not (p / "FAILED.txt").exists():
            continue
        rows.append((p.name, m, g, (p / "FAILED.txt").exists()))
    parts = [(f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>mscard gallery</title><style>{CSS}</style></head><body><h1>mscard: {len(rows)} subjects</h1>"
             "<div class='box'>Research software, not a medical device. Each row is one longitudinal pair; the dots are the calibration lights "
             "(new-lesion sensitivity, false new, false change on unchanged lesions, false resolved, SEL candidates on unchanged lesions, brain volume recovery) "
             "and the seventh dot is the real pair's own boundary check on brain volume (full mask against eroded mask).</div>"
             "<table><tr><th>Subject</th><th class='num'>Years</th><th class='num'>Baseline lesions</th><th class='num'>New</th><th class='num'>Enlarging</th>"
             "<th class='num'>Resolved</th><th class='num'>SEL cand.</th><th class='num'>Lesion vol %/yr</th><th class='num'>Brain vol %/yr</th>"
             "<th>Calibration</th><th class='num'>New found</th><th class='num'>Brain recovery</th></tr>")]
    for name, m, g, failed in rows:
        les = (m or {}).get("lesions", {})
        br = (m or {}).get("brain", {}).get("jacobian", {})
        L = (g or {}).get("lights", {})
        dots = "".join(_dot(L.get(k, "none")) for k in ("new_sensitivity", "false_new", "false_change_rate", "false_resolved_rate", "sel_candidates_per_lesion", "atrophy_recovery"))
        dots += _dot(br.get("boundary_light", "none"))
        n = (g or {}).get("new_lesions", {})
        parts.append(f"<tr><td><a href='{html.escape(name)}/report.html'>{html.escape(name)}</a>{' FAILED' if failed else ''}</td>"
                     f"<td class='num'>{_f((m or {}).get('interval_years'), 2)}</td><td class='num'>{_i(les.get('baseline_count'))}</td>"
                     f"<td class='num'>{_i(les.get('new'))}</td><td class='num'>{_i(les.get('enlarging'))}</td><td class='num'>{_i(les.get('resolved'))}</td>"
                     f"<td class='num'>{_i((m or {}).get('sel', {}).get('candidates'))}</td><td class='num'>{_f(les.get('volume_change_pct_per_year'), 1)}</td>"
                     f"<td class='num'>{_f(br.get('change_pct_per_year'), 2)}</td><td>{dots}</td>"
                     f"<td class='num'>{_i(n.get('reported'))}/{_i(n.get('injected'))}</td>"
                     f"<td class='num'>{_f((g or {}).get('brain', {}).get('jacobian_recovery'), 2)}</td></tr>")
    parts.append(f"</table><p class='small'>mscard {html.escape(__version__)}, generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}.</p></body></html>")
    p = out / "index.html"
    p.write_text("\n".join(parts))
    return p
