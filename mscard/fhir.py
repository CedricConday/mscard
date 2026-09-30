"""mscard.json (+ grade.json) -> FHIR R4B Bundle.

One DiagnosticReport (status ``preliminary``: research software) and one Observation per finding.
Each Observation carries the calibration light as ``interpretation`` (a local code system, since no
standard code expresses "recovered by the same pipeline on this scan") and the calibration evidence
as a note. Condition code: SNOMED CT 24700007 (multiple sclerosis). Finding codes are local, under
CODESYSTEM, because no LOINC code exists for them.
"""

from __future__ import annotations

import json
import uuid
from pathlib import Path

from fhir.resources.R4B.bundle import Bundle

CODESYSTEM = "https://github.com/CedricConday/mscard/fhir/CodeSystem/finding"
LIGHTS = "https://github.com/CedricConday/mscard/fhir/CodeSystem/calibration-light"
FINDINGS = [  # code, display, path into mscard.json, unit, light key
    ("new-lesions", "New T2 lesions", ("lesions", "new"), "{lesions}", "new_sensitivity"),
    ("enlarging-lesions", "Enlarging T2 lesions", ("lesions", "enlarging"), "{lesions}", "false_change_rate"),
    ("shrinking-lesions", "Shrinking T2 lesions", ("lesions", "shrinking"), "{lesions}", "false_change_rate"),
    ("resolved-lesions", "Resolved T2 lesions", ("lesions", "resolved"), "{lesions}", "false_resolved_rate"),
    ("lesion-volume-change", "T2 lesion volume change", ("lesions", "volume_change_pct_per_year"), "%/a", "false_change_rate"),
    ("sel-candidates", "Slowly expanding lesion candidates", ("sel", "candidates"), "{lesions}", "sel_candidates_per_lesion"),
    ("brain-volume-change", "Brain volume change (deformation field)", ("brain", "jacobian", "change_pct_per_year"), "%/a", "atrophy_recovery"),
]


def _get(d, path):
    for k in path:
        d = (d or {}).get(k)
    return d


def bundle(mscard: dict, grade: dict | None = None, patient_ref: str | None = None) -> dict:
    subj = {"reference": patient_ref or f"Patient/{mscard['subject']}"}
    lights = (grade or {}).get("lights", {})
    obs, refs = [], []
    for code, disp, path, unit, lk in FINDINGS:
        v = _get(mscard, path)
        if v is None:
            continue
        oid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"{mscard['subject']}/{mscard['follow_up']}/{code}"))
        light = lights.get(lk, "none")
        o = {"resourceType": "Observation", "id": oid, "status": "preliminary",
             "code": {"coding": [{"system": CODESYSTEM, "code": code, "display": disp}], "text": disp},
             "subject": subj, "valueQuantity": {"value": round(float(v), 3), "unit": unit, "system": "http://unitsofmeasure.org", "code": unit},
             "interpretation": [{"coding": [{"system": LIGHTS, "code": light}],
                                 "text": {"green": "read at face value", "amber": "read with the calibration",
                                          "red": "not trustworthy on this scan", "none": "not graded"}[light]}],
             "note": [{"text": f"mscard {mscard['mscard_version']}; interval {mscard['interval_years']:.2f} years "
                               f"({mscard['baseline']} to {mscard['follow_up']})"}]}
        obs.append(o)
        refs.append({"reference": f"Observation/{oid}"})
    rid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"{mscard['subject']}/{mscard['follow_up']}/report"))
    report = {"resourceType": "DiagnosticReport", "id": rid, "status": "preliminary",
              "code": {"text": "Longitudinal MS MRI report with self-calibration (mscard)"},
              "subject": subj, "result": refs,
              "conclusion": "Research software, not a medical device. Each result carries the recovery of injected known "
                            "change on this patient's own baseline as its interpretation.",
              "conclusionCode": [{"coding": [{"system": "http://snomed.info/sct", "code": "24700007",
                                              "display": "Multiple sclerosis"}]}]}
    b = {"resourceType": "Bundle", "type": "collection",
         "entry": [{"fullUrl": f"urn:uuid:{r['id']}", "resource": r} for r in [report, *obs]]}
    Bundle.model_validate(b)  # raises on any structural error
    return b


def from_dir(subject_dir: Path, patient_ref: str | None = None) -> dict:
    subject_dir = Path(subject_dir)
    m = json.loads((subject_dir / "measure" / "mscard.json").read_text())
    g = subject_dir / "calibration" / "grade.json"
    return bundle(m, json.loads(g.read_text()) if g.exists() else None, patient_ref)
