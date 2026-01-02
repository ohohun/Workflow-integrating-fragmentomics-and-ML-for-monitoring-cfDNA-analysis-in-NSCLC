import pandas as pd

configfile: "config/config.yaml"

# read sample sheet
df = pd.read_csv(config["samplesheet"], sep="\t")
req_cols = {"sample", "R1", "R2"}
missing = req_cols - set(df.columns)
if missing:
    raise ValueError(f"samples.tsv missing columns: {missing}")

SAMPLES = df["sample"].tolist()

rule all:
    input:
        expand("results/cnv/{s}.seg", s=SAMPLES),
        "results/reports/multiqc_report.html"

# 1) Trim adapters (fastp) —  .fastq or .fastq.gz
rule trim_fastp:
    input:
        r1=lambda w: df.set_index("sample").loc[w.s, "R1"],
        r2=lambda w: df.set_index("sample").loc[w.s, "R2"]
    output:
        r1="results/trim/{s}_R1.trimmed.fastq.gz",
        r2="results/trim/{s}_R2.trimmed.fastq.gz",
        html="results/trim/{s}_fastp.html",
        json="results/trim/{s}_fastp.json"
    threads: 8
    shell:
        r"""
        mkdir -p results/trim
        R2="{input.r2}"
        if [ -z "$R2" ] || [ "$R2" = "NA" ]; then
            fastp -i {input.r1} -o {output.r1} -h {output.html} -j {output.json} -w {threads}
            ln -sf {output.r1} {output.r2}
        else
            fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} -w {threads}
        fi
        """

# 2) FastQC หลัง Trim (เช็ค R1 เป็นหลัก)
rule fastqc:
    input:
        r1 = "results/trim/{s}_R1.trimmed.fastq.gz"
    output:
        html = "results/fastqc/{s}_R1.trimmed_fastqc.html",
        zip  = "results/fastqc/{s}_R1.trimmed_fastqc.zip"
    threads: 2
    shell:
        r"""
        mkdir -p results/fastqc
        fastqc {input.r1} -o results/fastqc
        """

# 3) Mapping ด้วย bwa-mem2 → SAM
rule map_bwa_mem2:
    input:
        r1="results/trim/{s}_R1.trimmed.fastq.gz",
        r2="results/trim/{s}_R2.trimmed.fastq.gz"
    output:
        temp("results/bam/{s}.sam")
    threads: config["threads_map"]
    params:
        ref=config["reference_fa"],
        extra=config.get("bwa_extra",""),
        rg=lambda w: "@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA".format(s=w.s)
    shell:
        r"""
        mkdir -p results/bam
        if [ "$(readlink -f {input.r1})" = "$(readlink -f {input.r2})" ]; then
            bwa-mem2 mem -t {threads} {params.extra} -R "{params.rg}" {params.ref} {input.r1} > {output}
        else
            bwa-mem2 mem -t {threads} {params.extra} -R "{params.rg}" {params.ref} {input.r1} {input.r2} > {output}
        fi
        """

# 4) Sort BAM
rule sort_bam:
    input: "results/bam/{s}.sam"
    output: temp("results/bam/{s}.sorted.bam")
    threads: config["threads_sort"]
    shell:
        r"""
        samtools sort -@ {threads} -O bam -o {output} {input}
        """

# 5) Mark duplicates (Picard)
rule mark_duplicates:
    input: "results/bam/{s}.sorted.bam"
    output:
        bam="results/bam/{s}.markdup.bam",
        metrics="results/bam/{s}.dup_metrics.txt"
    shell:
        r"""
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=SILENT
        """

# 6) Index BAM
rule index_bam:
    input: "results/bam/{s}.markdup.bam"
    output: "results/bam/{s}.markdup.bam.bai"
    shell:
        r"""
        samtools index {input}
        """

# 7) ichorCNA
rule make_wig:
    input:
        bam="results/bam/{s}.markdup.bam"
    output:
        wig="results/cnv/{s}.wig"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/cnv
        CHRS=$(bash -lc 'for i in $(seq 1 22); do printf "chr%s," "$i"; done; echo -n "chrX,chrY"')
        readCounter -w 1000000 -q 20 --chromosome "$CHRS" {input.bam} > {output.wig}
        """

rule ichorcna:
    input:
        wig="results/cnv/{s}.wig",
        gc="ref/hg38/gc_hg38_1000kb.wig",
        map="ref/hg38/map_hg38_1000kb.wig"
    output:
        "results/cnv/{s}.seg"
    params:
        outdir="results/cnv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        Rscript --vanilla bin/runIchorCNA.R \
          --id {wildcards.s} \
          --WIG {input.wig} \
          --gcWig {input.gc} \
          --mapWig {input.map} \
          --normal 'c(0.5,0.7,0.9)' \
          --ploidy 'c(2)' \
          --chrTrain 'c(paste0("chr",1:22),"chrX")' \
          --genomeStyle UCSC \
          --outDir {params.outdir}
        mv {params.outdir}/{wildcards.s}*.seg {output} 2>/dev/null || true
        """
      


# 8) MultiQC

rule multiqc:
    input:
        expand("results/fastqc/{s}_R1.trimmed_fastqc.html", s=SAMPLES),
        expand("results/trim/{s}_fastp.html", s=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        """
        mkdir -p results/multiqc
        multiqc results -o results/multiqc
        """

