import os, glob, re

FR = config["fragmentomics"]
OUT = FR["outdir"]

import csv

SAMPLES_TSV = "config/samples.tsv"

with open(SAMPLES_TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    SAMPLES = sorted({row["sample"] for row in reader})

def bam_clean(wc):
    return FR["bam_clean"].format(sample=wc.sample)

rule fragmentomics_all:
    input:
        expand(os.path.join(OUT, "{sample}", "fraglen", "fraglen.tsv"), sample=SAMPLES),
        expand(os.path.join(OUT, "{sample}", "endmotif", "endmotif_k4.tsv"), sample=SAMPLES),
        expand(os.path.join(OUT, "{sample}", "delfi", "delfi.tsv"), sample=SAMPLES),
        os.path.join(OUT, "fragmentomics_features.tsv")

rule make_autosome_chrom_sizes:
    input:
        fai=FR["fai"]
    output:
        chrom_sizes=FR["chrom_sizes"]
    shell:
        r"""
        mkdir -p $(dirname {output.chrom_sizes})
        cut -f1,2 {input.fai} \
          | awk -F'\t' -v OFS='\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {{print $1,$2}}' \
          > {output.chrom_sizes}
        """

rule make_delfi_bins_100kb:
    input:
        chrom_sizes=FR["chrom_sizes"]
    output:
        bins=FR["bins_100kb"]
    conda:
        "../envs/finaletoolkit.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bins})
        bedtools makewindows -g {input.chrom_sizes} -w 100000 > {output.bins}
        """

rule make_gap_bed:
    output:
        gap=FR["gap_bed"]
    shell:
        r"""
        mkdir -p $(dirname {output.gap})
        finaletoolkit gap-bed hg38 {output.gap}
        """

rule finale_fraglen_bins:
    input:
        bam=bam_clean
    output:
        tsv=os.path.join(OUT, "{sample}", "fraglen", "fraglen.tsv"),
        hist=os.path.join(OUT, "{sample}", "fraglen", "fraglen.hist.png")
    params:
        min_len=FR["fraglen"]["min"],
        max_len=FR["fraglen"]["max"],
        bin_size=FR["fraglen"]["bin_size"],
        short_fraction=FR["fraglen"]["short_fraction"],
        mapq=FR["mapq"]
    conda:
        "../envs/finaletoolkit.yaml"
    benchmark:
        os.path.join("benchmarks", "finale_fraglen_bins", "{sample}.tsv")
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        finaletoolkit frag-length-bins \
          -min {params.min_len} -max {params.max_len} \
          --bin-size {params.bin_size} \
          -stats -sf {params.short_fraction} \
          --histogram-path {output.hist} \
          -q {params.mapq} \
          -o {output.tsv} \
          {input.bam}
        """

rule finale_end_motifs:
    input:
        bam=bam_clean,
        twobit=FR["twobit"]
    output:
        tsv=os.path.join(OUT, "{sample}", "endmotif", "endmotif_k4.tsv")
    params:
        k=FR["endmotif"]["k"],
        min_len=FR["endmotif"]["min"],
        max_len=FR["endmotif"]["max"],
        mapq=FR["mapq"],
        workers=FR["workers"],
        both_strands=FR["endmotif"]["both_strands"]
    conda:
        "../envs/finaletoolkit.yaml"
    threads:
        FR["workers"]
    benchmark:
        os.path.join("benchmarks", "finale_end_motifs", "{sample}.tsv")
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        EXTRA=""
        if [ "{params.both_strands}" = "false" ]; then
          EXTRA="-B"
        fi

        finaletoolkit end-motifs \
          -k {params.k} -min {params.min_len} -max {params.max_len} \
          $EXTRA \
          -q {params.mapq} -w {params.workers} \
          -o {output.tsv} \
          {input.bam} {input.twobit}
        """

rule finale_delfi:
    input:
        bam=bam_clean,
        chrom_sizes=FR["chrom_sizes"],
        twobit=FR["twobit"],
        bins=FR["bins_100kb"],
        gap=FR["gap_bed"]          # ✅ เพิ่มบรรทัดนี้
    output:
        tsv=os.path.join(OUT, "{sample}", "delfi", "delfi.tsv")
    params:
        mapq=FR["mapq"],
        workers=FR["workers"],
        window_size=FR["delfi"]["window_size"],
        merge_bins=FR["delfi"]["merge_bins"],
        gc_correct=FR["delfi"]["gc_correct"]
    conda:
        "../envs/finaletoolkit.yaml"
    threads:
        FR["workers"]
    benchmark:
        os.path.join("benchmarks", "finale_delfi", "{sample}.tsv")
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})

        EXTRA=""
        if [ "{params.merge_bins}" = "false" ]; then
          EXTRA="$EXTRA -M"
        fi
        if [ "{params.gc_correct}" = "false" ]; then
          EXTRA="$EXTRA -G"
        fi

        finaletoolkit delfi \
          -g {input.gap} \
          -s {params.window_size} \
          -q {params.mapq} -w {params.workers} \
          $EXTRA \
          -o {output.tsv} \
          {input.bam} {input.chrom_sizes} {input.twobit} {input.bins}
        """

rule merge_fragmentomics_features:
    input:
        fraglen=expand(os.path.join(OUT, "{sample}", "fraglen", "fraglen.tsv"), sample=SAMPLES),
        endmotif=expand(os.path.join(OUT, "{sample}", "endmotif", "endmotif_k4.tsv"), sample=SAMPLES),
        delfi=expand(os.path.join(OUT, "{sample}", "delfi", "delfi.tsv"), sample=SAMPLES)
    output:
        tsv=os.path.join(OUT, "fragmentomics_features.tsv")
    conda:
        "../envs/finaletoolkit.yaml"
    script:
        "scripts/merge_fragmentomics_features.py"

