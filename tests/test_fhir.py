import pytest
from pydantic import ValidationError

from mscard.fhir import bundle

M = {"subject": "P1", "mscard_version": "0.1.0", "interval_years": 1.4, "baseline": "T1", "follow_up": "T3",
     "lesions": {"new": 4, "enlarging": 3, "shrinking": 4, "resolved": 7, "volume_change_pct_per_year": -18.0},
     "sel": {"candidates": 20}, "brain": {"jacobian": {"change_pct_per_year": 1.3}}}
G = {"lights": {"new_sensitivity": "none", "false_change_rate": "green", "sel_candidates_per_lesion": "red", "atrophy_recovery": "amber"}}


def test_bundle_valid_and_complete():
    b = bundle(M, G)
    kinds = [e["resource"]["resourceType"] for e in b["entry"]]
    assert kinds.count("Observation") == 7 and kinds[0] == "DiagnosticReport"
    sel = next(e["resource"] for e in b["entry"] if e["resource"].get("code", {}).get("coding", [{}])[0].get("code") == "sel-candidates")
    assert sel["interpretation"][0]["coding"][0]["code"] == "red"


def test_invalid_rejected():
    bad = dict(M, subject=None, lesions={"new": "x"})
    with pytest.raises((ValidationError, TypeError, ValueError, KeyError)):
        bundle(bad, G)
