"""Read-only integrity and numerical checks for the pre-deposit E9 package.

This script checks the supplied files. Cross-assay specimen identity was
confirmed from laboratory records by the study authors; it cannot be tested
from a de-identified CSV. This script does not recreate the historical P test.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from asca_revision import FIELDS, data_cube, decompose


ROOT = Path(__file__).resolve().parent


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)
    print(f"PASS | {message}")


def main() -> None:
    matrix = pd.read_csv(ROOT / "asca_input.csv", dtype={"strain": str})
    keys = ["strain", "treatment", "replicate"]
    require(matrix.shape == (60, 15), "60 design rows, 12 endpoint columns")
    require(set(FIELDS) == set(matrix.columns) - set(keys), "12 expected endpoint names")
    require(not matrix[keys].duplicated().any(), "unique design-cell replicate keys")
    require(matrix.groupby(["strain", "treatment"]).size().eq(6).all(), "six ordinals per design cell")
    require(np.isfinite(matrix[FIELDS].to_numpy(dtype=float)).all(), "no missing or non-finite endpoint values")

    cag = matrix.loc[(matrix.strain == "26695") & (matrix.treatment == "Blank NP")
                     & (matrix.replicate == 3), "Phospho_total_CagA"]
    require(len(cag) == 1 and np.isclose(cag.iloc[0], 1.248), "author-confirmed CagA record 1.248")

    batches = pd.read_csv(ROOT / "formulation_batches.csv")
    require(len(batches) == 15 and batches.groupby("formulation").size().eq(3).all(), "15 batch records")
    contrasts = pd.read_csv(ROOT / "adjusted_contrasts.csv")
    require(len(contrasts) == 541 and "Source" in contrasts.columns, "541 labelled comparisons")
    require(not contrasts["Source"].isna().any(), "all adjusted comparisons have provenance labels")
    require(len(pd.read_csv(ROOT / "reported_primary_comparisons.csv")) == 38,
            "38 separately transcribed reported primary comparisons")
    quality = pd.read_csv(ROOT / "exclusion_inventory.csv")
    require(len(quality) == 8 and not quality[["Exclusion_count", "Exclusion_reason"]].isna().any().any(),
            "exclusion reporting status is explicit for all eight assay categories")

    ct = pd.read_csv(ROOT / "reference_ct_lab_first_series.csv")
    require(len(ct) == 60 and not ct[["Strain", "Treatment", "Culture ordinal"]].duplicated().any(),
            "60 unique reference Ct culture-ordinal records")
    for gene in ["gyrB", "16S"]:
        cols = [f"{gene} Ct {i}" for i in (1, 2, 3)]
        require(np.isfinite(ct[cols].to_numpy(dtype=float)).all(), f"all {gene} Ct triplicates present")

    primers = pd.read_csv(ROOT / "primer_efficiency_status.csv")
    row = primers.loc[primers.Assay == "16S rRNA"]
    require(len(primers) == 8 and len(row) == 1, "eight recorded primer assays")
    require(np.isclose(row.Measured_efficiency_pct.iloc[0], 99.6), "16S measured efficiency: 99.6%")
    expected = {"cagA": 99.8, "vacA": 100.0, "ureA": 99.8, "flaA": 99.7,
                "babA": 99.5, "sabA": 99.6, "gyrB": 99.7, "16S rRNA": 99.6}
    require(all(np.isclose(primers.set_index("Assay").loc[g, "Measured_efficiency_pct"], val)
                for g, val in expected.items()), "all eight measured efficiencies match the supplied validation records")

    inventory = pd.read_csv(ROOT / "source_inventory.csv")
    for rec in inventory.to_dict("records"):
        filename = ROOT / "archive_csv" / rec["archive_path"]
        require(filename.is_file(), f"source present: {rec['archive_path']}")
        digest = hashlib.sha256(filename.read_bytes()).hexdigest()
        require(digest == rec["sha256"], f"source digest matches: {rec['archive_path']}")

    _, ss = decompose(data_cube(ROOT / "asca_input.csv"))
    fraction = ss / ss.sum() * 100
    target = np.array([52.482604576977415, 4.269889334839547,
                       2.9127699931471924, 40.33473609503584])
    require(np.allclose(fraction, target, atol=1e-8), "ASCA partition matches independent implementation")
    old_result = ROOT / "results" / "summary.json"
    if old_result.is_file():
        summary = json.loads(old_result.read_text(encoding="utf-8"))
        require(np.isclose(summary["sensitivity_P"]["interaction"], 0.5216478352164784),
                "bundled sensitivity output identifies P = 0.521648")
    rebuilt = ROOT / "rebuilt" / "asca_input.csv"
    if rebuilt.is_file():
        check = pd.read_csv(rebuilt, dtype={"strain": str})
        require(matrix.equals(check), "recreated matrix matches archived matrix exactly")

    print("\nNumerical/data-file checks passed. Matched specimen identity is based on the study's laboratory records, not inferable from de-identified CSV alone.")


if __name__ == "__main__":
    main()
