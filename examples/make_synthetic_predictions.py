"""Create fictitious test scores, never surrogate evidence for the manuscript."""

import argparse
import csv
import math
import random
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    rng = random.Random(20261029)
    rows = []
    for group in ("SYNTHETIC_TRAIN", "SYNTHETIC_VALIDATION"):
        for i in range(150):
            x1, x2 = rng.gauss(0, 1), rng.gauss(0, 1)
            sigmoid = lambda x: 1 / (1 + math.exp(-x))
            truth = sigmoid(-2 + .8 * x1 + .7 * x2)
            rows.append({"patient_id": f"SYN_{len(rows) + 1:04}", "cohort": group,
                "outcome": int(rng.random() < truth),
                "clinical_probability": sigmoid(-2 + .8 * x1),
                "radiomics_probability": sigmoid(-2 + .7 * x2),
                "fusion_probability": truth})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    print("Wrote 300 synthetic records; no real patient data.")


if __name__ == "__main__":
    main()
