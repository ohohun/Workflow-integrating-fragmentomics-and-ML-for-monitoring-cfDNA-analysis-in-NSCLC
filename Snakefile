import pandas as pd

configfile: "config/config.yaml"
df = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = df["sample"].tolist()

rule all:
    input:
        expand("results/cnv/{s}.seg", s=SAMPLES),
        "results/reports/multiqc_report.html"

# Trim adapters
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
        if [ "$R2" = "NA" ] || [ -z "$R2" ]; then
            fastp -i {input.r1} -o {output.r1} -h {output.html} -j {output.json} -w {threads}
        else
            fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} -w {threads}
        fi
        """

# QC หลัง Trim
rule fastqc:
    input:
        r1="results/trim/{s}_R1.trimmed.fastq.gz"
    output:
        "results/fastqc/{s}_R1_fastqc.html"
    shell:
        "mkdir -p results/fastqc; fastqc {input} -o results/fastqc"

# Mapping
rule map_bwa_mem2:
    input:
        r1="results/trim/{s}_R1.trimmed.fastq.gz",
        r2="results/trim/{s}_R2.trimmed.fastq.gz"
    output:
        temp("results/bam/{s}.sam")
    threads: config["threads_map"]
    params:
        ref=config["reference_fa"],
        extra=config["bwa_extra"],
        rg=lambda w: "@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA".format(s=w.s)
    shell:
        r"""
        mkdir -p results/bam
        if [ ! -f {input.r2} ]; then
            bwa-mem2 mem -t {threads} {params.extra} -R "{params.rg}" {params.ref} {input.r1} > {output}
        else
            bwa-mem2 mem -t {threads} {params.extra} -R "{params.rg}" {params.ref} {input.r1} {input.r2} > {output}
        fi
        """

# Sort BAM
rule sort_bam:
    input: "results/bam/{s}.sam"
    output: temp("results/bam/{s}.sorted.bam")
    threads: config["threads_sort"]
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input}"

# Mark duplicates
rule mark_duplicates:
    input: "results/bam/{s}.sorted.bam"
    output:
        bam="results/bam/{s}.markdup.bam",
        metrics="results/bam/{s}.dup_metrics.txt"
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=SILENT"

# Index BAM
rule index_bam:
    input: "results/bam/{s}.markdup.bam"
    output: "results/bam/{s}.markdup.bam.bai"
    shell:
        "samtools index {input}"

# ichorCNA CNV
rule ichorcna:
    input:
        bam="results/bam/{s}.markdup.bam",
        bai="results/bam/{s}.markdup.bam.bai"
    output: "results/cnv/{s}.seg"
    params:
        ref=config["reference_fa"],
        gc=config["gc_wig"],
        mapwig=config["map_wig"],
        binsize=config["binsize_cna"],
        ploidy=",".join([str(x) for x in config["ploidy"]]),
        normal=config["normal"],
        outdir="results/cnv"
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript -e 'library(ichorCNA);
            ichorCNA(bamFile="{input.bam}",
                     gcWig="{params.gc}",
                     mapWig="{params.mapwig}",
                     normal="{params.normal}",
                     ploidy=c({params.ploidy}),
                     refFasta="{params.ref}",
                     outDir="{params.outdir}",
                     id="{wildcards.s}",
                     binsize={params.binsize},
                     includeHOMD=FALSE,
                     maxCN=5,
                     chrTrain=c(paste0("chr",1:22),"chrX"),
                     writeSegmentTable=TRUE)'
        mv {params.outdir}/{wildcards.s}*.seg {output}
        """

# MultiQC Summary
rule multiqc:
    input:
        expand("results/fastqc/{s}_R1_fastqc.html", s=SAMPLES),
        expand("results/trim/{s}_fastp.html", s=SAMPLES)
    output: "results/reports/multiqc_report.html"
    shell:
        "multiqc -o results/reports ."
