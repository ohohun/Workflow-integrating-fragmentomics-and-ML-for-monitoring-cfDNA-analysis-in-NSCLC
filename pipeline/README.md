# cfDNA NSCLC Bioinformatics Pipeline

This folder contains the bioinformatics pipeline for processing low-coverage whole-genome sequencing (lcWGS) cell-free DNA (cfDNA) data for non-small cell lung cancer (NSCLC) monitoring.

The pipeline is designed to convert raw FASTQ sequencing files into quality control reports, BAM files, copy number alteration (CNA) features, tumor fraction features, and fragmentomics features that can be used for downstream visualization and machine learning analysis.

---

## Project Overview

**Project title:**  
Workflow integrating fragmentomics and machine learning for monitoring cell-free DNA analysis in non-small cell lung cancer

This pipeline supports cfDNA-based NSCLC monitoring by extracting features from lcWGS cfDNA data, including:

- Copy number alteration profiles
- Tumor fraction outputs from ichorCNA
- Fragment length distribution
- End-motif profiles
- MultiQC quality-control summary reports
- Merged fragmentomics feature table for downstream analysis

---

## Folder Structure

```text
pipeline/
├── Snakefile
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── envs/
│   └── *.yaml
├── rules/
│   ├── fragmentomics_finale.smk
│   └── other workflow rule files
├── scripts/
│   └── merge_fragmentomics_features.py
├── multiqc_config.yaml
├── ichor_targets.txt
└── README.md
```

### Folder Description

| Path | Description |
|---|---|
| `Snakefile` | Main Snakemake workflow file |
| `config/` | Configuration files and sample sheet |
| `config/config.yaml` | Main parameter and path configuration |
| `config/samples.tsv` | Input sample sheet containing FASTQ file paths |
| `envs/` | Conda environment files for workflow tools |
| `rules/` | Modular Snakemake rule files |
| `scripts/` | Python scripts for feature processing and merging |
| `multiqc_config.yaml` | MultiQC configuration file |
| `ichor_targets.txt` | Target/chromosome information for ichorCNA-related analysis |
| `README.md` | Pipeline documentation |

---

## Main Technologies

| Technology / Tool | Purpose |
|---|---|
| **Snakemake** | Workflow management and reproducible pipeline execution |
| **Python** | Feature processing, merging, and downstream data handling |
| **R** | Running ichorCNA and related CNA analysis |
| **fastp** | Adapter trimming and read quality filtering |
| **FastQC** | Sequencing read quality assessment |
| **MultiQC** | Aggregation of QC reports across samples |
| **BWA-MEM2** | Alignment of cfDNA sequencing reads to the human reference genome |
| **Samtools** | BAM file sorting, indexing, and manipulation |
| **Picard** | Marking duplicate reads |
| **HMMcopy readCounter** | Read counting in genomic bins for CNV analysis |
| **ichorCNA** | Copy number alteration analysis and tumor fraction estimation from cfDNA |
| **FinaleToolkit** | cfDNA fragmentomics analysis |
| **Bedtools** | Genomic interval and window processing |
| **Conda** | Environment and dependency management |

---

## Input Requirements

### 1. FASTQ Files

The pipeline accepts paired-end FASTQ files:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

Single-end data can also be used by setting `R2` as `NA` in the sample sheet.

FASTQ files are not included in this repository because sequencing files are usually large and may contain sensitive information.

---

### 2. Sample Sheet

Edit the sample sheet:

```text
config/samples.tsv
```

Example for paired-end data:

```tsv
sample	R1	R2
S01	data/fastq/S01_R1.fastq.gz	data/fastq/S01_R2.fastq.gz
S02	data/fastq/S02_R1.fastq.gz	data/fastq/S02_R2.fastq.gz
```

Example for single-end data:

```tsv
sample	R1	R2
S01	data/fastq/S01.fastq.gz	NA
```

Required columns:

| Column | Description |
|---|---|
| `sample` | Sample name used for output files |
| `R1` | Path to read 1 FASTQ file |
| `R2` | Path to read 2 FASTQ file, or `NA` for single-end data |

---

## Configuration

Before running the workflow, edit:

```text
config/config.yaml
```

Important paths to check:

```yaml
samplesheet: "config/samples.tsv"

reference_fa: "ref/hg38/hg38.fa"

fragmentomics:
  fasta: "ref/hg38/hg38.fa"
  fai: "ref/hg38/hg38.fa.fai"
  twobit: "ref/hg38/hg38.2bit"
  chrom_sizes: "ref/hg38/hg38.autosomes.chrom.sizes"
  bins_100kb: "ref/hg38/hg38.autosomes.100kb.bed"
  gap_bed: "ref/hg38/hg38.gap.bed"
```

Make sure all reference files exist before running the pipeline.

---

## Required Reference Files

The pipeline requires human genome reference files, preferably hg38:

```text
ref/hg38/hg38.fa
ref/hg38/hg38.fa.fai
ref/hg38/hg38.2bit
ref/hg38/gc_hg38_1000kb.wig
ref/hg38/map_hg38_1000kb.wig
```

Additional derived files may be generated or required for fragmentomics analysis:

```text
ref/hg38/hg38.autosomes.chrom.sizes
ref/hg38/hg38.autosomes.100kb.bed
ref/hg38/hg38.gap.bed
```

---

## How to Run Pipeline

### 1. Clone the Repository

```bash
git clone https://github.com/ohohun/Workflow-integrating-fragmentomics-and-ML-for-monitoring-cfDNA-analysis-in-NSCLC.git
cd Workflow-integrating-fragmentomics-and-ML-for-monitoring-cfDNA-analysis-in-NSCLC
```

---

### 2. Go to the Pipeline Folder

```bash
cd pipeline
```

---

### 3. Prepare FASTQ Files

Place or link FASTQ files in a local directory, for example:

```text
data/fastq/
├── S01_R1.fastq.gz
├── S01_R2.fastq.gz
├── S02_R1.fastq.gz
└── S02_R2.fastq.gz
```

---

### 4. Edit the Sample Sheet

Open and edit:

```bash
nano config/samples.tsv
```

Example:

```tsv
sample	R1	R2
S01	data/fastq/S01_R1.fastq.gz	data/fastq/S01_R2.fastq.gz
S02	data/fastq/S02_R1.fastq.gz	data/fastq/S02_R2.fastq.gz
```

---

### 5. Edit the Configuration File

Open and edit:

```bash
nano config/config.yaml
```

Check that these paths match your local environment:

```yaml
samplesheet: "config/samples.tsv"
reference_fa: "ref/hg38/hg38.fa"
```

Also check the `fragmentomics` reference paths in the same file.

---

### 6. Create and Activate Conda Environment

Create a Conda environment:

```bash
conda create -n cfdna_pipeline -c conda-forge -c bioconda \
  snakemake fastp fastqc multiqc bwa-mem2 samtools picard bedtools python r-base
```

Activate the environment:

```bash
conda activate cfdna_pipeline
```

Install FinaleToolkit:

```bash
pip install finaletoolkit
```

Make sure `ichorCNA`, `HMMcopy`, and `readCounter` are installed and accessible from the command line.

---

### 7. Check Required Tools

```bash
which snakemake
which fastp
which fastqc
which multiqc
which bwa-mem2
which samtools
which bedtools
which finaletoolkit
```

If any command is not found, install the missing tool before running the workflow.

---

### 8. Dry Run

Run a dry run to check whether the workflow is correctly configured:

```bash
snakemake -n
```

This command does not run the pipeline. It only checks the workflow, input files, output files, and rule dependencies.

---

### 9. Run the Full Pipeline

Run the workflow using 8 CPU cores:

```bash
snakemake --cores 8
```

Run with Snakemake-managed Conda environments:

```bash
snakemake --use-conda --cores 8
```

For a local computer with fewer resources, reduce the number of cores:

```bash
snakemake --use-conda --cores 4
```

---

### 10. Run Selected Outputs

Generate only the MultiQC report:

```bash
snakemake results/multiqc/multiqc_report.html --cores 8
```

Generate CNV segmentation output for one sample:

```bash
snakemake results/cnv/S01.seg --cores 8
```

Generate the merged fragmentomics feature table:

```bash
snakemake results/fragmentomics/fragmentomics_features.tsv --cores 8
```

---

## Pipeline Workflow

```text
FASTQ files
   ↓
Adapter trimming and read filtering
   ↓
Quality control
   ↓
Alignment to hg38 reference genome
   ↓
BAM sorting and duplicate
   ↓
BAM indexing
   ↓
CNV and tumor fraction analysis
   ↓
Fragmentomics analysis
   ↓
Feature merging
   ↓
Final outputs for ML and visualization
```

---

## Main Analysis Steps

### 1. Read Preprocessing

Tool: `fastp`

Main outputs:

```text
results/trim/{sample}_R1.trimmed.fastq.gz
results/trim/{sample}_R2.trimmed.fastq.gz
results/trim/{sample}_fastp.html
results/trim/{sample}_fastp.json
```

Purpose:

- Remove adapter sequences
- Filter low-quality reads
- Generate preprocessing reports

---

### 2. Quality Control

Tools: `FastQC`, `MultiQC`

Main outputs:

```text
results/fastqc/
results/multiqc/multiqc_report.html
```

Purpose:

- Evaluate sequencing read quality
- Review GC content
- Review duplication level
- Summarize QC across all samples

---

### 3. Genome Alignment

Tools: `BWA-MEM2`, `Samtools`, `Picard`

Main outputs:

```text
results/bam/{sample}.sorted.bam
results/bam/{sample}.markdup.bam
results/bam/{sample}.markdup.bam.bai
```

Purpose:

- Align cfDNA reads to the human reference genome
- Sort BAM files
- Mark duplicate reads
- Index final BAM files for downstream analysis

---

### 4. CNV and Tumor Fraction Analysis

Tools: `HMMcopy readCounter`, `ichorCNA`

Main outputs:

```text
results/cnv/{sample}.wig
results/cnv/{sample}.seg
```

Purpose:

- Count reads in genomic bins
- Correct GC and mappability bias
- Detect copy number alteration patterns
- Support tumor fraction interpretation from cfDNA

---

### 5. Fragmentomics Analysis

Tool: `FinaleToolkit`

Main analyses:

```text
Fragment length distribution
End-motif analysis
DELFI-style genome-wide fragmentation profile
```

Main outputs:

```text
results/fragmentomics/{sample}/fraglen/fraglen.tsv
results/fragmentomics/{sample}/fraglen/fraglen.hist.png
results/fragmentomics/{sample}/endmotif/endmotif_k4.tsv
results/fragmentomics/{sample}/delfi/delfi.tsv
```

Purpose:

- Characterize cfDNA fragment size patterns
- Analyze short and long fragment profiles
- Extract 4-mer end-motif signatures
- Generate genome-wide fragmentation profiles

---

### 6. Feature Merging

Tool: `Python`

Main output:

```text
results/fragmentomics/fragmentomics_features.tsv
```

Purpose:

- Merge fragment length, end-motif, and DELFI-style outputs
- Create a sample-level feature table
- Prepare features for machine learning and web-based visualization

---

## Final Outputs

After a successful run, the main outputs are:

```text
results/
├── trim/
│   ├── *_R1.trimmed.fastq.gz
│   ├── *_R2.trimmed.fastq.gz
│   ├── *_fastp.html
│   └── *_fastp.json
├── fastqc/
│   └── *_fastqc.html
├── bam/
│   ├── *.sorted.bam
│   ├── *.markdup.bam
│   └── *.markdup.bam.bai
├── cnv/
│   ├── *.wig
│   └── *.seg
├── multiqc/
│   └── multiqc_report.html
└── fragmentomics/
    ├── {sample}/
    │   ├── fraglen/
    │   ├── endmotif/
    │   └── delfi/
    └── fragmentomics_features.tsv
```

---

## Output Interpretation

| Output | Description |
|---|---|
| `multiqc_report.html` | Summary report of sequencing and preprocessing QC |
| `.markdup.bam` | Final aligned BAM file after duplicate marking |
| `.markdup.bam.bai` | BAM index file |
| `.wig` | Read count signal for ichorCNA |
| `.seg` | Copy number segmentation result |
| `fraglen.tsv` | Fragment length distribution |
| `fraglen.hist.png` | Fragment length histogram |
| `endmotif_k4.tsv` | 4-mer end-motif profile |
| `delfi.tsv` | Genome-wide fragmentation profile |
| `fragmentomics_features.tsv` | Merged feature table for ML and visualization |

---

## Common Errors and Solutions

### Missing FASTQ files

Check paths in:

```text
config/samples.tsv
```

Use:

```bash
ls data/fastq/
```

---

### Missing reference genome files

Check:

```bash
ls ref/hg38/
```

Required files may include:

```text
hg38.fa
hg38.fa.fai
hg38.2bit
gc_hg38_1000kb.wig
map_hg38_1000kb.wig
```

---

### Snakemake command not found

Install Snakemake:

```bash
conda install -c conda-forge -c bioconda snakemake
```

---

### FinaleToolkit path issue

If the workflow contains machine-specific paths, update:

```text
rules/fragmentomics_finale.smk
```

Example:

```python
FTK_BIN = "finaletoolkit"
BEDTOOLS = "bedtools"
SAMTOOLS = "samtools"
```

---

## Recommended `.gitignore`

Large sequencing files, reference genomes, and generated results should not be uploaded to GitHub.

Recommended `.gitignore` entries:

```gitignore
results/
logs/
ref/
*.fastq
*.fastq.gz
*.bam
*.bai
*.sam
*.wig
*.seg
*.html
*.pdf
.DS_Store
__pycache__/
```

---

## Example Complete Run

```bash
cd pipeline

conda activate cfdna_pipeline

snakemake -n

snakemake --use-conda --cores 8
```

---

## Notes

- This pipeline is intended for research use.
- The workflow is designed for low-coverage WGS cfDNA data.
- Reference genome paths must be updated according to the local environment.
- Large sequencing files such as FASTQ, BAM, and reference genome files should not be committed to GitHub.
- Sensitive patient level data should not be uploaded to a public repository.

---

## Contact

This pipeline was developed as part of a Health Data Science senior project on cfDNA-based NSCLC monitoring.
