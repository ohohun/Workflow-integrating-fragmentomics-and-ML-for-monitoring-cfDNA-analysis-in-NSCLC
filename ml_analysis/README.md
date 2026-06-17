# Machine Learning Analysis

This folder contains the downstream cfDNA NSCLC machine learning analysis module. It is designed to be used after the bioinformatics `pipeline/` has generated cfDNA-derived feature tables.

## Folder Structure

```text
ml_analysis/
├── cfDNA_NSCLC_ML_analysis.py
├── requirements.txt
├── README.md
└── data/
    └── README.md
```

## Main Script

`cfDNA_NSCLC_ML_analysis.py` is the main analysis script. It integrates:

- ichorCNA-derived features, such as tumor fraction, ploidy, and GC MAD
- fragment length features
- end-motif features
- clinical and treatment outcome data

The script covers data loading, cleaning, feature engineering, exploratory analysis, statistical testing, machine learning evaluation, feature interpretation, and result export.

## Required Input Files

Place the following local files in `ml_analysis/data/` before running the analysis:

```text
outcome.xlsx
clinicalData.xlsx
CRA_fraglen.csv
CRA_ichorCNA.csv
CRA_endmotif.csv
EGA_fraglen.csv
EGA_ichorCNA.csv
EGA_endmotif.csv
```

Files exported from your computer may have names such as `CRA_fraglen(18).csv`. Rename them to the clean names above before running.

## Installation

```bash
cd ml_analysis
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows, activate the environment with:

```bash
.venv\Scripts\activate
```

## Run Analysis

```bash
python cfDNA_NSCLC_ML_analysis.py
```

By default, the script reads input files from:

```text
data/
```

and writes generated outputs to:

```text
outputs/
```

You can override these paths with environment variables:

```bash
CFDNA_DATA_DIR=/path/to/data CFDNA_OUTPUT_DIR=/path/to/outputs python cfDNA_NSCLC_ML_analysis.py
```

## Data Privacy

Real clinical data, patient-level metadata, raw sequencing data, and generated outputs should not be committed to GitHub. This repository only includes the analysis code and documentation.
