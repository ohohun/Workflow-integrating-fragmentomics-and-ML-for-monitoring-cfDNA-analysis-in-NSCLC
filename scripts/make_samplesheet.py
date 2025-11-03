#!/usr/bin/env python3
import sys, csv, pathlib, re

ROOT = pathlib.Path("data")
GROUPS = ["cancer", "healthy"]

def find_fastq_pairs(folder: pathlib.Path):
    # รองรับ .fastq และ .fastq.gz
    fastqs = list(folder.glob("*.fastq")) + list(folder.glob("*.fastq.gz"))
    r1s, r2s = {}, {}
    for f in fastqs:
        name = f.name
        base = re.sub(r"(_R[12]|_[12])\.fastq(\.gz)?$", "", name)
        if re.search(r"_R1\.fastq(\.gz)?$|_1\.fastq(\.gz)?$", name):
            r1s[base] = f
        elif re.search(r"_R2\.fastq(\.gz)?$|_2\.fastq(\.gz)?$", name):
            r2s[base] = f
    pairs = []
    for base in sorted(set(r1s) & set(r2s)):
        pairs.append((r1s[base], r2s[base]))
    return pairs

rows = []
for g in GROUPS:
    gdir = ROOT / g
    if not gdir.exists():
        continue
    for egaf in sorted([p for p in gdir.iterdir() if p.is_dir()]):
        pairs = find_fastq_pairs(egaf)
        if not pairs:
            sys.stderr.write(f"[warn] no FASTQ pairs in {egaf}\n")
            continue
        for i, (R1, R2) in enumerate(pairs, 1):
            sample = egaf.name if len(pairs) == 1 else f"{egaf.name}__set{i}"
            rows.append({
                "sample": sample,
                "type": g,
                "R1": str(R1),
                "R2": str(R2),
                "library": "",
                "lane": "",
                "center": ""
            })

with open("samples.tsv", "w", newline="") as fo:
    w = csv.DictWriter(
        fo,
        fieldnames=["sample","type","R1","R2","library","lane","center"],
        delimiter="\t"
    )
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"[ok] wrote samples.tsv with {len(rows)} entries")
