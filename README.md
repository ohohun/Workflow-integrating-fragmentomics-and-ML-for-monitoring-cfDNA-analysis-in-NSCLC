# Workflow Integrating Fragmentomics and Machine Learning for Monitoring cfDNA Analysis in NSCLC

This repository contains an integrated workflow for cell-free DNA (cfDNA) analysis in non-small cell lung cancer (NSCLC). The project combines bioinformatics processing, machine learning analysis, and web-based visualization to support research on cfDNA-based cancer monitoring.

The overall aim of this project is to develop a non-invasive workflow for monitoring treatment response in NSCLC using cfDNA sequencing data. The workflow integrates genomic features, fragmentomic patterns, machine learning models, and clinical interpretation tools to support continuous disease tracking and precision-based research.

---

## Project Overview

Cell-free DNA has potential as a minimally invasive biomarker for cancer monitoring because it can be collected from blood samples and analyzed repeatedly over time. In this project, low-coverage whole-genome sequencing (lcWGS) cfDNA data are processed to extract molecular features related to NSCLC, including copy number alterations, tumor fraction, fragment length profiles, and end-motif patterns.

The extracted features are then used for machine learning analysis and visualization through a web application. The repository is organized into three main components:

1. **Bioinformatics pipeline**
2. **Machine learning analysis**
3. **Web application for result visualization**

---

## Repository Structure

```text
Workflow-integrating-fragmentomics-and-ML-for-monitoring-cfDNA-analysis-in-NSCLC/
├── pipeline/
│   └── Bioinformatics workflow for processing cfDNA sequencing data
│
├── ml_analysis/
│   └── Machine learning analysis using extracted cfDNA features
│
├── web_app/
│   └── Web-based application for visualizing cfDNA analysis results
│
├── README.md
└── .gitignore
```

---

## Main Components

### 1. Bioinformatics Pipeline

The `pipeline/` folder contains the bioinformatics workflow used to process cfDNA sequencing data. This part focuses on transforming raw or preprocessed sequencing data into analysis-ready genomic and fragmentomic features.

The pipeline is designed to support key cfDNA analysis steps such as:

- sequencing data quality control
- read alignment and preprocessing
- copy number alteration analysis
- tumor fraction estimation
- fragment length feature extraction
- end-motif feature preparation
- generation of output files for downstream analysis

This component provides the foundation for extracting meaningful cfDNA-derived features that can be used in the machine learning and visualization steps.

---

### 2. Machine Learning Analysis

The `ml_analysis/` folder contains the machine learning analysis module. This part uses cfDNA-derived features to explore patterns associated with NSCLC and treatment response.

Main objectives of this component include:

- preparing feature tables for machine learning
- integrating cfDNA features from multiple sources
- performing feature selection
- training and evaluating machine learning models
- supporting classification or prediction tasks related to NSCLC monitoring
- generating model outputs for interpretation and visualization

Example feature groups used in the analysis include:

- copy number alteration features
- tumor fraction-related features
- fragment length distribution features
- short-to-long fragment ratios
- end-motif or motif-derived features
- longitudinal changes in cfDNA features

The main script in this folder is:

```text
ml_analysis/cfDNA_NSCLC_ML_analysis.py
```

---

### 3. Web Application

The `web_app/` folder contains a web-based application for visualizing cfDNA analysis results. This component is intended to make the output of the bioinformatics and machine learning workflows easier to interpret.

The web application can support visualization at both cohort and patient levels, including:

- patient-level cfDNA profiles
- longitudinal monitoring results
- tumor fraction trends
- genomic feature summaries
- machine learning prediction outputs
- interactive result exploration

This component helps connect the computational workflow with user-friendly result interpretation.

---

## Workflow Summary

The overall workflow of this project can be summarized as follows:

```text
cfDNA sequencing data
        ↓
Bioinformatics pipeline
        ↓
Feature extraction
        ↓
Machine learning analysis
        ↓
Prediction and interpretation
        ↓
Web-based visualization
```

---

## Main Technologies

This project uses a combination of bioinformatics, data science, and web development tools.

### Bioinformatics

- Python
- Snakemake
- BWA
- samtools
- ichorCNA
- MultiQC
- cfDNA fragmentomics-related analysis tools

### Machine Learning and Data Analysis

- Python
- pandas
- NumPy
- scikit-learn
- scipy
- matplotlib
- seaborn

### Web Application

- Python
- Flask or Streamlit-based workflow
- HTML
- CSS
- JavaScript
- data visualization libraries

---

## How to Run

### 1. Clone the Repository

```bash
git clone https://github.com/ohohun/Workflow-integrating-fragmentomics-and-ML-for-monitoring-cfDNA-analysis-in-NSCLC.git
cd Workflow-integrating-fragmentomics-and-ML-for-monitoring-cfDNA-analysis-in-NSCLC
```

---

### 2. Run the Bioinformatics Pipeline

Go to the pipeline folder:

```bash
cd pipeline
```

Install the required environment or dependencies based on the pipeline configuration.

Then run the workflow:

```bash
snakemake --cores 4
```

The pipeline will process cfDNA sequencing data and generate output files for downstream feature analysis.

---

### 3. Run Machine Learning Analysis

Go to the machine learning folder:

```bash
cd ../ml_analysis
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Run the analysis script:

```bash
python cfDNA_NSCLC_ML_analysis.py
```

Input data should be prepared and placed in the appropriate data folder before running the analysis.

---

### 4. Run the Web Application

Go to the web application folder:

```bash
cd ../web_app
```

Install the required dependencies according to the web application setup.

Then run the application. For example:

```bash
python app.py
```

or, if the application uses Streamlit:

```bash
streamlit run app.py
```

After running the application, open the local URL shown in the terminal to view the dashboard.

---

## Expected Outputs

This project can generate several types of outputs, including:

- cfDNA quality control reports
- copy number alteration profiles
- tumor fraction estimation results
- fragment length features
- end-motif features
- machine learning model performance results
- prediction outputs
- patient-level visualization dashboards
- cohort-level summary dashboards

---

## Data Availability

Patient-level sequencing and clinical data are not included in this repository due to privacy and data sharing restrictions.

Users should prepare their own input data and place the required files in the appropriate folders before running the pipeline or machine learning analysis.

---

## Project Purpose

This repository was developed as a research-support workflow for cfDNA-based monitoring in NSCLC. It is intended for academic and research use, especially for exploring how fragmentomics, machine learning, and visualization tools can be integrated to support non-invasive cancer monitoring.

---

## Notes

- The workflow is designed for research purposes only.
- The results should not be used as a standalone clinical decision-making tool.
- External validation and larger cohort testing are required before clinical application.
- Users should check each folder's README file for more specific instructions.

---

## Repository Modules

| Folder | Description |
|---|---|
| `pipeline/` | Bioinformatics workflow for processing cfDNA sequencing data and extracting genomic features |
| `ml_analysis/` | Machine learning analysis using cfDNA-derived features |
| `web_app/` | Web application for visualization and interpretation of cfDNA results |

---

## Author

This project was developed as part of a cfDNA-based NSCLC monitoring research workflow integrating bioinformatics, machine learning, and web-based visualization.
