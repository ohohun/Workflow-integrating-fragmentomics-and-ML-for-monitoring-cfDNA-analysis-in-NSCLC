import os
import pandas as pd
import numpy as np

def safe_read(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    return pd.read_csv(path, sep="\t")

def entropy(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return np.nan
    return float(-(p * np.log2(p)).sum())

def main():
    out_path = snakemake.output[0]
    samples = list(snakemake.params.samples)
    outdir = snakemake.params.outdir

    rows = []
    for s in samples:
        row = {"sample": s}

        frag_path = os.path.join(outdir, s, "fraglen", "fraglen.tsv")
        end_path  = os.path.join(outdir, s, "endmotif", "endmotif_k4.tsv")
        delfi_path = os.path.join(outdir, s, "delfi", "delfi.tsv")

        frag = safe_read(frag_path)
        if frag is not None:
            # เก็บ summary แบบปลอดภัย (ไม่ผูกชื่อคอลัมน์ตายตัว)
            row["fraglen_rows"] = frag.shape[0]
            for c in frag.columns:
                cl = c.lower()
                if cl in ("short_fraction","short_ratio","short_frag_ratio"):
                    row["fraglen_short_fraction"] = float(frag[c].iloc[0])

        endm = safe_read(end_path)
        if endm is not None:
            row["endmotif_rows"] = endm.shape[0]
            freq_col = None
            for c in endm.columns:
                if c.lower() in ("freq","frequency","fraction"):
                    freq_col = c
                    break
            if freq_col is not None:
                p = endm[freq_col].astype(float).values
                row["endmotif_entropy"] = entropy(p)
                row["endmotif_top1"] = float(np.nanmax(p))

        delfi = safe_read(delfi_path)
        if delfi is not None:
            num = delfi.select_dtypes(include="number")
            row["delfi_rows"] = delfi.shape[0]
            row["delfi_num_cols"] = num.shape[1]
            if num.shape[1] > 0:
                row["delfi_mean_of_means"] = float(num.mean().mean())

        rows.append(row)

    df = pd.DataFrame(rows).sort_values("sample")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
