"""Every threshold mscard uses, in one place, with its origin.

The traffic-light bands below are mscard's own working bands for reading a
calibration. They are not a clinical standard and no standard exists yet for
"how good must a longitudinal MS pipeline be"; they are stated so the reader
can disagree with them. Green means the pipeline recovered the injected truth
well enough that the corresponding number in the report can be read at face
value; amber means read it with the calibration next to it; red means the
number is not trustworthy on this scan with this pipeline.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SyntheticParams:
    volume_factor: float = 0.98  # whole-brain volume factor over the synthetic year, -2 %
    falloff_mm: float = 12.0  # bidsgate default: the contraction fades over this distance outside the brain
    n_new: int = 12  # new lesions injected, bidsgate's size cycle 30/80/200/600/1500 mm3
    exclusion_mm: float = 6.0  # new lesions are placed at least this far from every baseline lesion
    rigid_rot_deg: float = 1.0  # small head repositioning between visits (lesiontrack.synth defaults)
    rigid_trans_vox: float = 1.0
    intensity_scale: tuple = (0.9, 1.1)  # per-image gain drawn uniformly, models scanner drift
    noise_frac: float = 0.02  # Gaussian noise, fraction of the 99th percentile, inside the head only
    dt_years: float = 1.0
    seed: int = 0


@dataclass(frozen=True)
class Bands:
    """(green_at_least_or_at_most, amber_...) per grade; direction is in the name."""

    new_sensitivity_min: tuple = (0.90, 0.70)
    false_new_max: tuple = (1, 3)
    atrophy_recovery_range: tuple = ((0.8, 1.2), (0.5, 1.5))
    false_change_rate_max: tuple = (0.05, 0.15)  # baseline lesions called enlarging or shrinking when unchanged
    false_resolved_rate_max: tuple = (0.0, 0.05)
    sel_candidates_per_lesion_max: tuple = (0.2, 0.5)
    seg_sensitivity_min: tuple = (0.80, 0.60)  # only when mscard ran the segmenter itself
    seg_false_positive_max: tuple = (3, 8)


@dataclass(frozen=True)
class Params:
    synthetic: SyntheticParams = field(default_factory=SyntheticParams)
    bands: Bands = field(default_factory=Bands)
    threads: int = 4
    engine: str = "greedy"


LST_AI_IMAGE = "jqmcginnis/lst-ai:v2.0.0rc1-cpu"
