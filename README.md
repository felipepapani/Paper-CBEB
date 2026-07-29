# RQA of the ECG for Ventricular Fibrillation — SDDB

Reproducibility code for the CBEB 2026 manuscript:

> **A Short-Scale Transition into Ventricular Fibrillation: Amplitude and Heart-Rate Variability Recurrence Contrast**
> F. G. Papani, H. H. B. Nunes, D. C. Soriano, F. J. R. Sales.

This repository contains the analysis pipeline that applies Recurrence
Quantification Analysis (RQA) to electrocardiographic recordings from the
Sudden Cardiac Death Holter Database (SDDB), plus the reviewer-requested
control analyses.

## What is here

- `TCC_reproducibility.ipynb` — the full pipeline, from raw ECG to the figures
  and effect sizes reported in the paper.
- `R4_contraste_mesmos_7.py` — recomputes the amplitude-domain paired Cohen's *d*
  restricted to the seven HRV-arm patients (control for the domain contrast).
- `R5_verificar_ritmo.py` — reads the rhythm annotations in the final 30 s before
  the fibrillation onset (checks for ventricular tachycardia).
- `requirements.txt` — Python dependencies.

## What is NOT here

The **raw SDDB signals are not redistributed**. The database is publicly
available from PhysioNet under the Open Data Commons Attribution License (ODC-BY):

- https://physionet.org/content/sddb/1.0.0/

Download it from PhysioNet and point the notebook's `DATA_DIR` / Drive paths to
your local copy.

## Pipeline overview

1. Load results already saved to disk (Part 0)
2. Data acquisition — SDDB (Parts 1-2)
3. Core functions: NaN-robust filter, RQA (15 metrics), percentile binarization (Part 3)
4. Filter-selection test (Part 4)
5. Main RQA extraction over thresholds p05-p50 (Part 5)
6. Threshold and feature comparison (Part 6)
7. Predictive test — far (-60 to -30 min) vs near (-30 to 0 min) (Part 7)
8. Per-patient paired statistics (Part 8)
9. HRV / Takens recurrence arm (Part 9)
10. Reviewer analyses — R-4, R-5, gate coverage, bootstrap CIs (Part 10)
11. Figure regeneration with English labels (Part 11)

## Method summary

- Recurrence threshold: percentile of the positive pairwise distances per window
  (swept p05-p50; p40 adopted).
- Amplitude arm: embedding dimension m = 1, Euclidean distance, l_min = v_min = 2,
  5-second windows.
- HRV arm: NN series from audited sinus beats; Takens embedding with tau from the
  first minimum of average mutual information and m from the false-nearest-neighbors
  criterion; 256-beat windows, 32-beat step; 80% NN-density quality gate.
- Effect size: paired Cohen's d_z (far vs near), one value per patient per phase.

## Environment

```bash
pip install -r requirements.txt
```

Python 3.10+. The notebook was developed on Google Colab.

## Citation

If you use this code, please cite the CBEB 2026 paper and this repository
(the Zenodo DOI badge will appear here once the first release is archived).

## License

Code: MIT (see `LICENSE`). Data: not included; SDDB is ODC-BY via PhysioNet.
