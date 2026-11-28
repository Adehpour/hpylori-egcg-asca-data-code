# EGCG–chitosan/TPP nanoparticles and *H. pylori*: analysis data and code

Version 1.1.0. This repository contains the data and independently implemented statistical scripts supporting the two-isolate *H. pylori* EGCG nanoparticle study. The v1.1.0 update adds the 35-sheet final supplementary-data workbook and harmonizes the documentation with the revised manuscript; the underlying canonical CSV inputs and ASCA code remain unchanged.

## Licenses and citation

The research datasets, CSV exports, spreadsheet, derived results, figures and
documentation are released under **CC BY 4.0** (see `LICENSE-DATA.md`). The
Python source code is released under the **MIT License** (see `LICENSE`).
Please credit the study authors and cite the **v1.1.0** release using
`CITATION.cff` and its version-specific DOI,
**[10.5281/zenodo.24622353](https://doi.org/10.5281/zenodo.24622353)**.
The earlier **v1.0.0** version has DOI
**[10.5281/zenodo.24595291](https://doi.org/10.5281/zenodo.24595291)**.
For references to the complete version series, the concept DOI is
**[10.5281/zenodo.24595290](https://doi.org/10.5281/zenodo.24595290)**.
The versioned archive is also available via the
[GitHub v1.1.0 release](https://github.com/Adehpour/hpylori-egcg-asca-data-code/releases/tag/v1.1.0).

## Figure numbering in the revised manuscript

The final supplementary material contains one supplementary figure:
**Figure S1, paired 24-hour time-kill reductions** (called Figure S3 in an
earlier draft). The previously drafted reference-gene Ct and dual-reference
normalization plots (formerly Figures S1 and S2) are not included as figures
in the final supplementary material. Their underlying measurements and
normalization analyses remain available in `reference_ct_lab_first_series.csv`,
`Supplementary_Data_Final.xlsx`, and the other source tables
in this repository. Their removal as figures does not remove any data.

The plots in `results/Figure_7.png` and `results/Figure_7.pdf` are
independently generated **ASCA figures corresponding to manuscript Figure 8**;
their filenames are internal analysis outputs, not final figure numbers.
Likewise, figure labels within the original filenames under `archive_csv/`
identify historical assay exports and are not a numbering scheme for the
revised manuscript.

This figure-numbering note does not alter the canonical assay observations,
analysis code, or immutable v1.0.0 release archived under the DOI above.

## Dataset and provenance

The ASCA matrix (`asca_input.csv`) contains **60 biological samples, 12 measured endpoints, and three design columns**: strain, treatment and replicate number. The design has two isolates, five groups and six samples per group. **Laboratory records confirm that the 12 endpoint measurements in each row belong to the same physical bacterial culture.** Replicate numbers 1–6 are within-group sample labels; individual laboratory specimen identifiers are not necessary to interpret the anonymized matrix.

All 60 rows are complete and uniquely identified by the three design columns. No missing endpoint value has been imputed. `data_dictionary.csv` maps each variable to its unit, ASCA transformation and canonical source file. The 25 CSVs under `archive_csv/` are the authoritative input exports used by the reconstruction. Obsolete conflicting and redundant copies are excluded. In particular, the correct phospho-CagA/total-CagA value for strain 26695, blank nanoparticles, sample 3, is **1.248**.

| File | Contents |
| --- | --- |
| `asca_input.csv` | 60 matched biological samples × 12 numeric endpoints, with strain, treatment and replicate labels. |
| `data_dictionary.csv` | Variable definitions, input scales, transformations and source exports. |
| `archive_csv/` and `source_inventory.csv` | Canonical assay/formulation exports and SHA-256 checksums. |
| `reference_ct_lab_first_series.csv` | 60 culture records of gyrB and 16S rRNA technical Ct triplicates. |
| `primer_efficiency_status.csv` | All eight validated primer pairs and their **measured** amplification efficiencies. |
| `formulation_batches.csv` | Three physicochemical batch records for each of five nanoparticle formulations (15 total). |
| `adjusted_contrasts.csv` | 541 adjusted comparisons, with calculation family, data scale and source class. |
| `reported_primary_comparisons.csv` | 38 comparisons transcribed from the original manuscript, recorded separately from reanalyses. |
| `exclusion_inventory.csv` | Assay inclusion criteria and counts of analyzed records; not a historical laboratory incident ledger. |
| `time_kill.csv`, `mic_mbc.csv` | Available source-derived time-kill and dilution endpoint observations. |
| `results/` | ASCA variance partition, loadings, sample scores and independently generated figure. |
| `Supplementary_Data_Final.xlsx` | The 35-sheet supplementary workbook accompanying the final manuscript, including Tables S1–S10, adjusted-comparison reporting, formulation and microbiology records, primer validation, ASCA inputs and sensitivity analyses. This is the preferred workbook for article-specific tables and comparison labels. |
| `Supplementary_E9_Data_Code_Validated_v1.xlsx` | The original eight-sheet v1.0.0 reproducibility workbook, retained to preserve the earlier released package. The canonical machine-readable inputs remain the CSV files. |

The ASCA matrix and all 60 reference-gene Ct culture records agree numerically
between the two workbooks. The 35-sheet workbook expands the supplementary
tables and clarifies the provenance and statistical families of comparisons;
it does not replace or silently rewrite the original canonical assay exports.

The 541 contrasts are **not** represented as 541 original statistical-model printouts. They comprise 344 reconstructions from source CSVs, 180 exploratory reanalyses, 8 sensitivity comparisons using strain-normalized source values and 9 explicitly approximate comparisons reconstructed from published means/SDs. Keep the `Source`, `Method`, `Family` and `P_bound` columns when reusing this table. Analyses of normalized biofilm values are distinct from originally reported tests on raw absorbance; no raw OD values have been manufactured.

Blank numeric effect cells in `reported_primary_comparisons.csv` indicate that the corresponding effect estimate was not reported in the manuscript source; the reported P value and its comparison family remain recorded. `exclusion_inventory.csv` uses **Not reported** where the source material does not specify a retrospective exclusion count or reason. These labels do not mean that a missing value is zero.

## RT-qPCR validation

Amplification efficiencies, measured for the eight assays, are: *cagA* 99.8%, *vacA* 100.0%, *ureA* 99.8%, *flaA* 99.7%, *babA* 99.5%, *sabA* 99.6%, *gyrB* 99.7%, and 16S rRNA 99.6%. Primer sequences and literature sources are in `primer_efficiency_status.csv`. These results are within the reported 90–110% acceptance interval. They are actual reported measurements, not values inferred from the acceptance range.

## Computational methods

`asca_revision.py` implements a balanced ASCA partition of the 60 × 12 matrix. The six target-minus-*gyrB* ΔCt variables are negated, IL-8 is log10-transformed, and the other five outcomes retain their recorded scales. Each column is then mean-centered and autoscaled (sample SD). Biofilm is expressed as **percentage of untreated within the same strain** before the joint autoscaling. Its loading cannot be interpreted as absolute between-strain biomass.

The resulting total standardized variance partition is: **treatment 52.4826%, isolate 4.2699%, treatment × isolate 2.9128%, residual 40.3347%**. PC1 accounts for **94.5575% of treatment-associated variation**. The included Python code conducts 10,000 reproducible permutations (seed 271828) with a plus-one correction; the interaction sensitivity test produces **P = 0.521648**. Treatment and isolate each yield **P = 0.00009999**. A nonsignificant interaction with only two isolates does not demonstrate equivalence. The older interaction P = 0.481 came from a different historical analysis; the executable, documented implementation and its P = 0.521648 are the reproducible analysis provided with this release. Any manuscript reporting the released implementation should use the matching test specification and P value.

The included Python code is an **independent, documented implementation**. It is not falsely identified as the unavailable historical script. Source data, mathematical transformation, permutation strategy and random seed are sufficient to rerun the archived analysis.

## Reproduce

Use Python 3.12 with the pinned packages in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python verify_package.py
python asca_revision.py --input asca_input.csv --out results --seed 271828 --permutations 10000
python rebuild_tables.py
python verify_package.py
```

The first verification checks row counts, sample keys, source-file checksums, all eight efficiencies and variance partition. The ASCA script regenerates `results/summary.json`, `results/loadings.csv`, `results/scores.csv` and the ASCA figure; the table script regenerates the matrix and the computable comparison tables under `rebuilt/`. Generated output does not alter source exports. The paired sample identity is supported by the authors' laboratory confirmation and is not reverse-inferred from anonymized CSV values.

The original excluded/repeated-sample event log and original raw biofilm absorbance readings are not distributed. `exclusion_inventory.csv` specifies the available analyzed records and the stated criteria, without inventing a retrospective exclusion count or underlying measurements. The distributed data and code remain usable without those additional lab records.
