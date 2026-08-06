import os
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    from diffprivlib.mechanisms import Laplace
    _HAS_DIFFPRIVLIB = True
except Exception:
    _HAS_DIFFPRIVLIB = False

DATA_DIR = "data"
MISSING = {"--", "", "00", "XX", "NAN", "NONE"}
_ALLELE_CACHE = {}


def _laplace(value, epsilon, sensitivity=1.0):
    """Laplace mechanism. Falls back to numpy if diffprivlib is unavailable."""
    if _HAS_DIFFPRIVLIB:
        return Laplace(epsilon=epsilon, sensitivity=sensitivity).randomise(value)
    return value + np.random.default_rng().laplace(0.0, sensitivity / epsilon)


def build_allele_map(rsids):
    """Collect alleles observed at each rsID across the local genome cohort."""
    key = frozenset(rsids)
    if key in _ALLELE_CACHE:
        return _ALLELE_CACHE[key]

    alleles = defaultdict(set)
    files = [f for f in os.listdir(DATA_DIR)
             if f.startswith("genome_") and f.endswith(".txt")]
    for fname in files:
        with open(os.path.join(DATA_DIR, fname), "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 4 or p[0] not in rsids:
                    continue
                for a in p[3].strip().upper():
                    if a in "ACGT":
                        alleles[p[0]].add(a)

    _ALLELE_CACHE[key] = alleles
    return alleles


def apply_dp_noise(df: pd.DataFrame, epsilon: float = 1.0) -> pd.DataFrame:
    """
    Apply Laplace noise to LD-neighbor genotypes.

    Genotypes are encoded as allele dosage (0, 0.5, 1) relative to the two
    alleles observed at that locus, perturbed, then decoded back to a genotype
    built from those same alleles. This keeps output biologically valid: no
    allele is introduced that was never observed at the position.

    Rows are skipped when the genotype is missing, malformed, or when fewer
    than two alleles are observed at the locus (nothing to perturb toward).
    """
    df = df.copy()
    mask = df["ld_neighbor"] == True
    targets = df.loc[mask, "rsid"].tolist()
    if not targets:
        print("[DP] No LD neighbors tagged, nothing to perturb.")
        return df

    allele_map = build_allele_map(set(targets))

    perturbed = changed = skipped = 0

    for idx in df.index[mask]:
        rsid = df.at[idx, "rsid"]
        gt = str(df.at[idx, "genotype"]).strip().upper()

        if gt in MISSING or len(gt) != 2 or not all(a in "ACGT" for a in gt):
            skipped += 1
            continue

        observed = set(allele_map.get(rsid, set())) | set(gt)
        if len(observed) < 2:
            skipped += 1
            continue

        a1, a2 = sorted(observed)[:2]
        if not set(gt) <= {a1, a2}:
            skipped += 1
            continue

        dose = sum(1 for a in gt if a == a2) / 2.0
        noisy = _laplace(dose, epsilon)
        new_dose = int(round(min(max(noisy * 2, 0), 2)))
        new_gt = [a1 + a1, a1 + a2, a2 + a2][new_dose]

        df.at[idx, "genotype"] = new_gt
        perturbed += 1
        if new_gt != gt:
            changed += 1

    backend = "diffprivlib" if _HAS_DIFFPRIVLIB else "numpy fallback"
    print(f"[DP] epsilon={epsilon} ({backend}): perturbed {perturbed}, "
          f"changed {changed}, skipped {skipped}")
    return df