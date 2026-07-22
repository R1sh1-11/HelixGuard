import pandas as pd
from diffprivlib.mechanisms import Laplace

GENOTYPE_ENCODING = {
    "AA": 0.0, "AC": 0.5, "AG": 0.5, "AT": 0.5,
    "CC": 1.0, "CG": 0.5, "CT": 0.5,
    "GG": 1.0, "GT": 0.5, "TT": 1.0,
}


def _decode_genotype(noisy_val: float, original: str) -> str:
    best = min(GENOTYPE_ENCODING, key=lambda g: abs(GENOTYPE_ENCODING[g] - noisy_val))
    return best


def apply_dp_noise(df: pd.DataFrame, epsilon: float = 1.0) -> pd.DataFrame:
    df = df.copy()
    mask = df["ld_neighbor"] == True

    mech = Laplace(epsilon=epsilon, sensitivity=1.0)

    def noisy_genotype(row):
        encoded = GENOTYPE_ENCODING.get(row["genotype"], 0.5)
        noisy = mech.randomise(encoded)
        return _decode_genotype(noisy, row["genotype"])

    df.loc[mask, "genotype"] = df[mask].apply(noisy_genotype, axis=1)
    return df