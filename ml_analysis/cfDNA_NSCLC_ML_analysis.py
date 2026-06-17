# cfDNA NSCLC Machine Learning Analysis
#
# Converted from the original Jupyter notebook for GitHub-ready execution.
# Run from this folder:
#     python cfDNA_NSCLC_ML_analysis.py


# %% [cell 0]
# # cfDNA NSCLC Analysis Pipeline
#
# This notebook contains the analysis workflow used in the NSCLC cfDNA monitoring project. The code is organized from data loading and exploratory analysis through feature engineering, machine learning evaluation, model interpretation, and progression-free survival analysis.
#
# The main input path is controlled by `BASE`. If the data files are stored in the `data/` folder, set `BASE = "./data/"` in the load-data cell.
#
# ## 0 · Imports & Global Config

# %% [cell 1]
# %matplotlib inline  # Jupyter magic removed for .py execution
import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.model_selection import (StratifiedKFold, LeaveOneOut, cross_val_predict, learning_curve)
from sklearn.metrics import (
    roc_auc_score, roc_curve, accuracy_score, f1_score,
    confusion_matrix, ConfusionMatrixDisplay, recall_score,
    precision_score, brier_score_loss
)
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
from sklearn.feature_selection import mutual_info_classif
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from IPython.display import display
import warnings, os
from pathlib import Path
import lightgbm as lgb
import shap
warnings.filterwarnings('ignore')

# ── Global Reproducibility ────────────────────────────────────────────────────
SEED = 42

def set_global_seed(seed=SEED):
    """Set all random seeds for full reproducibility."""
    import random, os
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    try:
        import lightgbm as lgb
        lgb.LGBMClassifier()  # just to ensure lgb is loaded
    except Exception:
        pass

set_global_seed(SEED)

# ── Plot Style Configuration ─────────────────────────────────────────────────
plt.rcParams.update({
    'figure.dpi': 100,
    'figure.facecolor': 'white',
    'axes.facecolor': '#fafafa',
    'axes.grid': True,
    'grid.alpha': 0.25,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'font.family': 'DejaVu Sans',
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'legend.framealpha': 0.85,
})
sns.set_theme(style='whitegrid', palette='muted')
sns.set_context('notebook')

OUT = './'
os.makedirs(OUT, exist_ok=True)

# ── Colour palettes ──────────────────────────────────────────────
PAL = {
    'PR':      '#4C78A8',
    'SD':      '#F58518',
    'PD':      '#54A24B',
    'Healthy': '#888888',
    'NSCLC':   '#E45756',
}
C = {
    'PR':      '#2ecc71',
    'SD':      '#f39c12',
    'PD':      '#e74c3c',
    'NSCLC':   '#e74c3c',
    'Healthy': '#888888',   # FIX: เพิ่ม key ที่หายไป
}
TP_LABEL    = {1: 'Baseline\n(T0)', 2: '(W3-4)', 3: '(W12)'}
GROUP_ORDER = ['PR', 'SD', 'PD']

# ── Helper functions ─────────────────────────────────────────────
def pstar(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return f'ns (p={p:.2f})'

def draw_ref_band(ax, healthy_vals, color='#888888', alpha=0.15, label='Healthy range'):
    q1, q3 = np.percentile(healthy_vals, [25, 75])
    ax.axhspan(q1, q3, color=color, alpha=alpha, zorder=0)
    ax.axhline(np.median(healthy_vals), color=color, lw=1.2, ls='--',
               alpha=0.6, label=label)

def annotate_kw(ax, groups_vals, x_positions, y_ref, fontsize=10):
    vals = [v for v in groups_vals if len(v) > 1]
    if len(vals) < 2: return None
    _, p = kruskal(*vals)
    mid = np.mean(x_positions)
    ax.text(mid, y_ref, pstar(p), ha='center', va='bottom',
            fontsize=fontsize, fontweight='bold',
            color='#c0392b' if p < 0.05 else '#666666')
    return p


def add_panel_labels(axes_list, labels=None, fontsize=13, fontweight='bold',
                     x=-0.06, y=1.08):
    """Add (a)(b)(c)... labels to each subplot axis."""
    import string
    if labels is None:
        labels = list(string.ascii_lowercase)
    if hasattr(axes_list, 'flat'):
        axes_list = list(axes_list.flat)
    elif not isinstance(axes_list, (list, tuple)):
        axes_list = [axes_list]
    for i, ax in enumerate(axes_list):
        if i >= len(labels):
            break
        ax.text(x, y, f'({labels[i]})',
                transform=ax.transAxes,
                fontsize=fontsize, fontweight=fontweight,
                va='top', ha='right', color='black')

print('✓ Libraries loaded & config set')

# %% [cell 2]
# ## 1 · Load Data

# %% [cell 3]
# ── Project Paths ─────────────────────────────────────────────────────────────
# Run this notebook from the ml_analysis/ folder.
# Input files are expected inside ml_analysis/data/ by default.
# You can override paths with environment variables:
#   CFDNA_DATA_DIR=/path/to/data CFDNA_OUTPUT_DIR=/path/to/outputs
BASE_DIR = Path(os.getenv("CFDNA_DATA_DIR", "data"))
OUT_DIR  = Path(os.getenv("CFDNA_OUTPUT_DIR", "outputs"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Keep string variables for compatibility with the original notebook code.
BASE = str(BASE_DIR) + os.sep
OUT  = str(OUT_DIR) + os.sep

required_files = [
    "outcome.xlsx",
    "clinicalData.xlsx",
    "CRA_fraglen.csv",
    "CRA_ichorCNA.csv",
    "CRA_endmotif.csv",
    "EGA_fraglen.csv",
    "EGA_ichorCNA.csv",
    "EGA_endmotif.csv",
]

missing_files = [f for f in required_files if not (BASE_DIR / f).exists()]
if missing_files:
    raise FileNotFoundError(
        "Missing required input file(s) in ./data/: " + ", ".join(missing_files)
    )

# ── Load Data ────────────────────────────────────────────────────────────────
outcome   = pd.read_excel(BASE_DIR / "outcome.xlsx")
clinical  = pd.read_excel(BASE_DIR / "clinicalData.xlsx")
CRA_frag  = pd.read_csv(BASE_DIR / "CRA_fraglen.csv")
CRA_ichor = pd.read_csv(BASE_DIR / "CRA_ichorCNA.csv")
CRA_motif = pd.read_csv(BASE_DIR / "CRA_endmotif.csv")
EGA_frag  = pd.read_csv(BASE_DIR / "EGA_fraglen.csv")
EGA_ichor = pd.read_csv(BASE_DIR / "EGA_ichorCNA.csv")
EGA_motif = pd.read_csv(BASE_DIR / "EGA_endmotif.csv")

for nm, df in [('outcome', outcome), ('clinical', clinical), ('CRA_frag', CRA_frag),
               ('CRA_ichor', CRA_ichor), ('CRA_motif', CRA_motif),
               ('EGA_frag', EGA_frag), ('EGA_ichor', EGA_ichor),
               ('EGA_motif', EGA_motif)]:
    print(f'  {nm:12s}: {df.shape}')

# Add group labels
for df in [CRA_frag, CRA_ichor, CRA_motif, EGA_frag, EGA_ichor, EGA_motif]:
    df['label'] = df['disease_status'].map({'lung cancer': 'NSCLC', 'healthy': 'Healthy'})

all_frag  = pd.concat([CRA_frag,  EGA_frag],  ignore_index=True)
all_ichor = pd.concat([CRA_ichor, EGA_ichor], ignore_index=True)
all_motif = pd.concat([CRA_motif, EGA_motif], ignore_index=True)

# Feature column groups
ICHOR_FEATS = ['tumor_fraction', 'ploidy', 'gc_mad']
FRAG_FEATS  = ['frac_100_150', 'frac_151_180', 'frac_181_220', 'frac_221_300',
               'short_to_long_100_150_over_151_220', 'skewness',
               'kurtosis_excess', 'mean', 'median', 'short_fraction_s150']
MOTIF_COLS  = [c for c in CRA_motif.columns if len(c) == 4 and c.isalpha()]

print(f'\nOutcome: {outcome["CT"].value_counts().to_dict()}')
print(f'End-motif columns: {len(MOTIF_COLS)}')
print(f'Matching cancer PIDs in outcome: '
      f'{len([p for p in CRA_ichor[CRA_ichor["disease_status"]=="lung cancer"]["pid"].unique() if p in outcome["PID"].values])}')

# %% [cell 4]
# ---
# # Part 1 · EDA: What Does the cfDNA Data Tell Us?
#
# | Chapter | Question | Output |
# |---------|----------|--------|
# | 1 | ข้อมูลมีอะไรบ้าง? | Dataset Overview |
# | 2 | Healthy reference ค่าปกติ | CRA vs EGA Comparison |
# | 3 | Cancer vs Healthy ต่างกันแค่ไหน? | Baseline Discrimination |
# | 4 | Baseline — PR/SD/PD ต่างกันไหม? | Longitudinal Profiles |
# | 5 | Longitudinal trajectory | Trajectory by Response |
# | 6 | Fragmentomics & End-motif signal | Fragment + Motif Analysis |
# | 7 | Summary: feature ไหนสำคัญสุด? | Discriminability Heatmap |
#
# > 💡 **Interpretation Guide:** IQR band = Healthy reference range. Color coding: 🔵 PR, 🟠 SD, 🔴 PD

# %% [cell 5]
# ---
# ## Chapter 1 · Dataset Overview
#
# **Objective:** Understand sample composition, cohort distribution, and class balance.
#
# - **Panel A:** Total samples (NSCLC vs Healthy) across CRA + EGA cohorts  
# - **Panel B:** Sample distribution by timepoint  
# - **Panel C:** CRA cancer samples stratified by CT response (PR/SD/PD)

# %% [cell 6]
from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
def add_stacked_labels_general(
    ax, x, lower, upper,
    lower_name='CRA', upper_name='EGA',
    lower_inside_color='white', upper_inside_color='white',
    outside_color='#555555',
    small_thresh=18,
    total_offset=6,
    seg_offset=2.0,
    lower_fs=9, upper_fs=9, total_fs=10.5
):
    """
    Generic stacked-bar labels:
    - place inside segment if tall enough
    - otherwise place just above the segment
    """
    total = lower + upper

    # lower segment
    if lower > 0:
        if lower >= small_thresh:
            ax.text(
                x, lower / 2, f'{lower_name}\n{lower}',
                ha='center', va='center',
                fontsize=lower_fs, fontweight='bold',
                color=lower_inside_color, clip_on=False
            )
        else:
            ax.text(
                x, lower + seg_offset, f'{lower_name} {lower}',
                ha='center', va='bottom',
                fontsize=max(lower_fs - 0.5, 7.5),
                fontweight='bold', color=outside_color, clip_on=False
            )

    # upper segment
    if upper > 0:
        if upper >= small_thresh:
            ax.text(
                x, lower + upper / 2, f'{upper_name}\n{upper}',
                ha='center', va='center',
                fontsize=upper_fs, fontweight='bold',
                color=upper_inside_color, clip_on=False
            )
        else:
            ax.text(
                x, total + seg_offset, f'{upper_name} {upper}',
                ha='center', va='bottom',
                fontsize=max(upper_fs - 0.5, 7.5),
                fontweight='bold', color=outside_color, clip_on=False
            )

    # total label
    ax.text(
        x, total + total_offset, f'n={total}',
        ha='center', va='bottom',
        fontsize=total_fs, fontweight='bold',
        color='#333333', clip_on=False
    )


def add_simple_bar_labels(ax, bars, color, y_offset=0.45, fontsize=8.5):
    """Add numeric labels above bars."""
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + y_offset,
                f'{int(h)}',
                ha='center', va='bottom',
                fontsize=fontsize, fontweight='bold',
                color=color, clip_on=False
            )


# =========================================================
# Figure setup
# =========================================================
fig, axes = plt.subplots(1, 3, figsize=(18.8, 6.9))
fig.suptitle(
    'Dataset Overview: Sample Composition and Distribution',
    fontsize=13, fontweight='bold', y=0.98
)

# =========================================================
# Panel A: Sample Composition (CRA + EGA combined)
# =========================================================
ax = axes[0]
x_a = np.arange(2)
groups = ['NSCLC\n(CRA+EGA)', 'Healthy\n(CRA+EGA)']
cra_vals = [105, 10]
ega_vals = [42, 182]

ax.bar(
    x_a, cra_vals, width=0.52,
    color=['#E45756', '#888888'],
    alpha=0.90, edgecolor='white', linewidth=1.5, zorder=3
)
ax.bar(
    x_a, ega_vals, width=0.52, bottom=cra_vals,
    color=['#F4A261', '#BBBBBB'],
    alpha=0.90, edgecolor='white', linewidth=1.5, zorder=3
)

for i, (c, e) in enumerate(zip(cra_vals, ega_vals)):
    add_stacked_labels_general(
        ax, i, c, e,
        lower_name='CRA', upper_name='EGA',
        small_thresh=18,
        total_offset=5.5,
        seg_offset=2.0,
        lower_fs=9, upper_fs=9, total_fs=10.5
    )

ax.set_xticks(x_a)
ax.set_xticklabels(groups, fontsize=10)
ax.set_title(
    'A  |  Sample Composition: NSCLC versus Healthy\n(CRA + EGA Combined)',
    fontweight='bold', fontsize=10.5, pad=12
)
ax.set_ylabel('Samples (n)')
ax.set_ylim(0, 235)
ax.grid(axis='y', alpha=0.28, zorder=0)
ax.spines[['top', 'right']].set_visible(False)

ax.legend(
    handles=[
        Patch(facecolor='#E45756', alpha=0.90, label='CRA Cancer (105)'),
        Patch(facecolor='#F4A261', alpha=0.90, label='EGA Cancer (42)'),
        Patch(facecolor='#888888', alpha=0.90, label='CRA Healthy (10)'),
        Patch(facecolor='#BBBBBB', alpha=0.90, label='EGA Healthy (182)'),
    ],
    fontsize=7.2,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.01),
    framealpha=0.75,
    ncol=2,
    handlelength=1.2,
    columnspacing=1.0
)

ax.text(
    0.5, -0.18,
    'Total 339 samples · Cancer 147 (43%) · Healthy 192 (57%)',
    transform=ax.transAxes, ha='center',
    fontsize=8.5, color='#555555', style='italic'
)

# =========================================================
# Panel B: Sample distribution by timepoint
#   >>> custom-fixed to avoid EGA / n-label overlap
# =========================================================
ax = axes[1]

cra_nsclc_tp = CRA_ichor[CRA_ichor['disease_status'] == 'lung cancer'].groupby('timepoint').size()
ega_nsclc_tp = EGA_ichor[EGA_ichor['disease_status'] == 'lung cancer'].groupby('timepoint').size()
cra_hlt_tp = CRA_ichor[CRA_ichor['disease_status'] == 'healthy'].groupby('timepoint').size()
ega_hlt_tp = EGA_ichor[EGA_ichor['disease_status'] == 'healthy'].groupby('timepoint').size()

tps = [1, 2, 3]

# spaced out a bit more to reduce crowding
x_nsclc = np.array([1.0, 2.25, 3.50])
x_hlt = np.array([5.40])
divider_x = 4.35
bw_b = 0.68

cra_n = [cra_nsclc_tp.get(t, 0) for t in tps]
ega_n = [ega_nsclc_tp.get(t, 0) for t in tps]

# NSCLC bars
ax.bar(
    x_nsclc, cra_n, width=bw_b,
    color='#E45756', alpha=0.90,
    edgecolor='white', linewidth=1.5,
    label='NSCLC · CRA', zorder=3
)
ax.bar(
    x_nsclc, ega_n, width=bw_b, bottom=cra_n,
    color='#F4A261', alpha=0.90,
    edgecolor='white', linewidth=1.5,
    label='NSCLC · EGA', zorder=3
)

# Healthy reference bar
cra_h1 = cra_hlt_tp.get(1, 0)
ega_h1 = ega_hlt_tp.get(1, 0)

ax.bar(
    x_hlt, [cra_h1], width=bw_b,
    color='#888888', alpha=0.90,
    edgecolor='white', linewidth=1.5,
    label='Healthy · CRA', zorder=3
)
ax.bar(
    x_hlt, [ega_h1], width=bw_b, bottom=[cra_h1],
    color='#BBBBBB', alpha=0.90,
    edgecolor='white', linewidth=1.5,
    label='Healthy · EGA', zorder=3
)

# ---- custom labels for NSCLC bars (main fix) ----
for xi, c, e in zip(x_nsclc, cra_n, ega_n):
    # CRA segment
    if c >= 12:
        ax.text(
            xi, c / 2, f'CRA\n{c}',
            ha='center', va='center',
            fontsize=8.0, fontweight='bold',
            color='white', clip_on=False
        )
    else:
        ax.text(
            xi, c + 1.5, f'CRA {c}',
            ha='center', va='bottom',
            fontsize=7.5, fontweight='bold',
            color='#555555', clip_on=False
        )

    # EGA segment: always inside orange bar for these bars
    if e > 0:
        ax.text(
            xi, c + e / 2, f'EGA\n{e}',
            ha='center', va='center',
            fontsize=7.8, fontweight='bold',
            color='white', clip_on=False
        )

    # total label: lifted high enough so it won't collide with EGA label
    ax.text(
        xi, c + e + 7.0, f'n={c+e}',
        ha='center', va='bottom',
        fontsize=9.4, fontweight='bold',
        color='#333333', clip_on=False
    )

# ---- custom labels for Healthy bar ----
# Lower gray segment is small -> keep readable without colliding
if cra_h1 >= 12:
    ax.text(
        x_hlt[0], cra_h1 / 2, f'CRA\n{cra_h1}',
        ha='center', va='center',
        fontsize=8.0, fontweight='bold',
        color='white', clip_on=False
    )
else:
    ax.text(
        x_hlt[0], max(cra_h1 / 2, 6), f'CRA {cra_h1}',
        ha='center', va='center',
        fontsize=7.4, fontweight='bold',
        color='#4d4d4d', clip_on=False
    )

ax.text(
    x_hlt[0], cra_h1 + ega_h1 / 2, f'EGA\n{ega_h1}',
    ha='center', va='center',
    fontsize=8.0, fontweight='bold',
    color='white', clip_on=False
)

ax.text(
    x_hlt[0], cra_h1 + ega_h1 + 7.0, f'n={cra_h1+ega_h1}',
    ha='center', va='bottom',
    fontsize=9.4, fontweight='bold',
    color='#333333', clip_on=False
)

ax.set_xticks(list(x_nsclc) + list(x_hlt))
ax.set_xticklabels([TP_LABEL[t] for t in tps] + ['Healthy\n(TP1 only)'], fontsize=9)

ax.axvline(x=divider_x, color='#cccccc', lw=1.2, ls='--', zorder=1)
ax.text(
    divider_x + 0.05, 150, 'Healthy\nreference',
    fontsize=7.2, color='#888888',
    va='top', ha='left', style='italic',
    bbox=dict(facecolor='white', edgecolor='none', alpha=0.72, pad=1.4)
)

ax.set_title(
    'B  |  Sample Distribution by Timepoint\n(CRA + EGA Combined)',
    fontweight='bold', fontsize=10.5, pad=12
)
ax.set_ylabel('Samples (n)')
ax.set_ylim(0, 245)
ax.set_xlim(0.35, 6.15)
ax.grid(axis='y', alpha=0.28, zorder=0)
ax.spines[['top', 'right']].set_visible(False)

ax.legend(
    handles=[
        Patch(facecolor='#E45756', alpha=0.90, label='NSCLC · CRA'),
        Patch(facecolor='#F4A261', alpha=0.90, label='NSCLC · EGA'),
        Patch(facecolor='#888888', alpha=0.90, label='Healthy · CRA'),
        Patch(facecolor='#BBBBBB', alpha=0.90, label='Healthy · EGA'),
    ],
    fontsize=7.0,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.01),
    framealpha=0.75,
    ncol=2,
    handlelength=1.2,
    columnspacing=1.0
)

ax.text(
    0.5, -0.18,
    'Healthy available at TP1 only · NSCLC followed across 3 timepoints',
    transform=ax.transAxes, ha='center',
    fontsize=8.5, color='#555555', style='italic'
)

# =========================================================
# Panel C: CRA cancer samples by timepoint and response group
# =========================================================
ax = axes[2]

can_ichor_merged = (
    CRA_ichor[CRA_ichor['disease_status'] == 'lung cancer']
    .merge(outcome[['PID', 'CT']], left_on='pid', right_on='PID', how='left')
)

ct_groups = ['PR', 'SD', 'PD']
bw_c = 0.18
x_tp = np.array([1, 2, 3])
offsets_c = {'PR': -1.5 * bw_c, 'SD': -0.5 * bw_c, 'PD': 0.5 * bw_c}

for grp in ct_groups:
    vals = [
        can_ichor_merged[
            (can_ichor_merged['CT'] == grp) &
            (can_ichor_merged['timepoint'] == tp)
        ].shape[0]
        for tp in [1, 2, 3]
    ]

    bars_c = ax.bar(
        x_tp + offsets_c[grp], vals, width=bw_c,
        color=PAL[grp], alpha=0.88,
        edgecolor='white', linewidth=1.2,
        label=grp, zorder=3
    )

    add_simple_bar_labels(ax, bars_c, color=PAL[grp], y_offset=0.45, fontsize=8.5)

noct_data = [
    can_ichor_merged[
        (can_ichor_merged['CT'].isna()) &
        (can_ichor_merged['timepoint'] == tp)
    ].shape[0]
    for tp in [1, 2, 3]
]

bars_nc = ax.bar(
    x_tp + 1.5 * bw_c, noct_data, width=bw_c,
    color='#CCCCCC', alpha=0.80,
    edgecolor='#999999', linewidth=1.0,
    hatch='///', label='No CT', zorder=3
)

add_simple_bar_labels(ax, bars_nc, color='#999999', y_offset=0.45, fontsize=8.5)

for tp_i, tp in enumerate([1, 2, 3]):
    tot = (
        sum(
            can_ichor_merged[
                (can_ichor_merged['CT'] == g) &
                (can_ichor_merged['timepoint'] == tp)
            ].shape[0]
            for g in ct_groups
        )
        + noct_data[tp_i]
    )
    ax.text(
        tp, 24.0, f'Σ={tot}',
        ha='center', va='bottom',
        fontsize=8.5, color='#444444',
        fontweight='bold', clip_on=False
    )

ax.set_xticks([1, 2, 3])
ax.set_xticklabels([TP_LABEL[t] for t in [1, 2, 3]], fontsize=9)
ax.set_title(
    'C  |  CRA Cancer Samples by Timepoint and Response Group\n(PR / SD / PD / No CT Data)',
    fontweight='bold', fontsize=10.5, pad=12
)
ax.set_ylabel('Samples (n)')
ax.set_ylim(0, 27)
ax.grid(axis='y', alpha=0.28, zorder=0)
ax.spines[['top', 'right']].set_visible(False)

ax.legend(
    fontsize=8.2,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.01),
    framealpha=0.75,
    ncol=4,
    handlelength=1.2,
    columnspacing=1.0
)

ax.text(
    0.5, -0.18,
    '33 patients with CT outcome (PR+SD+PD) · 4 patients without CT',
    transform=ax.transAxes, ha='center',
    fontsize=8.5, color='#555555', style='italic'
)

# =========================================================
# Layout / Save / Show
# =========================================================
add_panel_labels(axes)
plt.tight_layout(rect=[0.02, 0.12, 0.98, 0.93], w_pad=2.8)
plt.savefig(f'{OUT}/eda_ch1_dataset_overview.png', dpi=150, bbox_inches='tight')
plt.show()

print('Key observations:')
print('  • A · Total 339 samples: NSCLC 147 (CRA 105 + EGA 42), Healthy 192 (CRA 10 + EGA 182)')
print('  • B · Healthy มีเฉพาะ TP1 (n=192); NSCLC ติดตาม 3 timepoints (51→49→47 samples)')
print('  • C · CRA cancer 37 patients: 33 มี CT outcome (PR=18, SD=9, PD=6), 4 ไม่มี CT')
print('  • Class imbalance PR:SD:PD = 18:9:6 → ต้องใช้ class_weight หรือ balanced metrics')

# %% [cell 7]
# ---
# ## Chapter 2 · Healthy Reference
#
# **Objective:** Establish normal cfDNA baseline values and check cohort consistency (CRA vs EGA).
#
# - Boxplots with jitter for `tumor_fraction`, `ploidy`, `gc_mad`
# - Mann–Whitney U test for CRA vs EGA difference
# - ⚠️ Cohort effect present if p < 0.05

# %% [cell 8]
hlt_cra = CRA_ichor[CRA_ichor['disease_status']=='healthy'].copy()
hlt_ega = EGA_ichor[EGA_ichor['disease_status']=='healthy'].copy()
hlt_cra['cohort'] = 'CRA Healthy'
hlt_ega['cohort'] = 'EGA Healthy'
hlt_all = pd.concat([hlt_cra, hlt_ega], ignore_index=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Healthy Reference: Comparison of CRA and EGA',
             fontsize=13, fontweight='bold', y=1.03)

for ax, feat in zip(axes, ICHOR_FEATS):
    data_c = hlt_cra[feat].dropna().values
    data_e = hlt_ega[feat].dropna().values
    _, p_mw = mannwhitneyu(data_c, data_e, alternative='two-sided')

    bp = ax.boxplot([data_c, data_e], patch_artist=True,
                    medianprops=dict(color='black', lw=2.2),
                    flierprops=dict(marker='o', markersize=4, alpha=0.5))
    bp['boxes'][0].set_facecolor('#888888'); bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor('#AAAAAA'); bp['boxes'][1].set_alpha(0.6)

    np.random.seed(42)
    for xi, vals in enumerate([data_c, data_e], 1):
        jx = xi + (np.random.rand(len(vals)) - 0.5) * 0.25
        ax.scatter(jx, vals, s=18, alpha=0.55, color='#555555',
                   edgecolors='none', zorder=4)

    ax.set_xticklabels([f'CRA\n(n={len(data_c)})', f'EGA\n(n={len(data_e)})'], fontsize=11)
    coh_str = pstar(p_mw) if p_mw < 0.05 else '(ns)'
    ax.set_title(f'{feat}\nMann–Whitney U, p = {p_mw:.3f}  {coh_str}', fontweight='bold', fontsize=11)
    ax.set_ylabel(feat, fontsize=10); ax.grid(axis='y', alpha=0.3)
    ax.spines[['top','right']].set_visible(False)
    ax.text(0.5, -0.22,
            f'CRA: {np.median(data_c):.4f} (IQR {np.percentile(data_c,25):.4f}–{np.percentile(data_c,75):.4f})\n'
            f'EGA: {np.median(data_e):.4f} (IQR {np.percentile(data_e,25):.4f}–{np.percentile(data_e,75):.4f})',
            transform=ax.transAxes, ha='center', fontsize=8, color='#555555')

add_panel_labels(axes)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/eda_ch2_healthy_reference.png', dpi=150, bbox_inches='tight')
plt.show()

print('Key observations:')
for feat in ICHOR_FEATS:
    dc = hlt_cra[feat].dropna(); de = hlt_ega[feat].dropna()
    _, p = mannwhitneyu(dc, de, alternative='two-sided')
    coh = 'cohort effect ⚠' if p < 0.05 else 'consistent ✓'
    print(f'  {feat:20s}: CRA median={dc.median():.4f}  EGA median={de.median():.4f}  ({coh})')

# %% [cell 9]
# ---
# ## Chapter 3 · Cancer vs Healthy
#
# **Objective:** Quantify cfDNA signal differences between NSCLC and healthy controls at baseline (TP1).
#
# - Gray band = Healthy IQR reference  
# - Optimal tumor_fraction threshold → used downstream as `OPT_TF`

# %% [cell 10]
tp1_ichor = all_ichor[all_ichor['timepoint']==1].copy()

fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
fig.suptitle('cfDNA Profiles: NSCLC versus Healthy Controls at Baseline (Timepoint1)',
             fontsize=14, fontweight='bold', y=1.01)

for ax, feat in zip(axes, ICHOR_FEATS):
    g_h = tp1_ichor[tp1_ichor['label']=='Healthy'][feat].dropna().values
    g_c = tp1_ichor[tp1_ichor['label']=='NSCLC'][feat].dropna().values
    _, p = mannwhitneyu(g_h, g_c, alternative='two-sided')

    draw_ref_band(ax, g_h, color='#888888', alpha=0.18,
                  label=f'Healthy IQR (n={len(g_h)})')

    bp = ax.boxplot([g_h, g_c], positions=[1, 2], patch_artist=True,
                    widths=0.45,
                    medianprops=dict(color='black', lw=2.2),
                    flierprops=dict(marker='o', ms=4, alpha=0.4))
    bp['boxes'][0].set_facecolor(PAL['Healthy']); bp['boxes'][0].set_alpha(0.55)
    bp['boxes'][1].set_facecolor(PAL['NSCLC']);   bp['boxes'][1].set_alpha(0.55)

    np.random.seed(7)
    for xi, vals, col in [(1, g_h, PAL['Healthy']), (2, g_c, PAL['NSCLC'])]:
        jx = xi + (np.random.rand(len(vals))-0.5)*0.3
        ax.scatter(jx, vals, s=16, alpha=0.5, color=col, edgecolors='none', zorder=4)

    ax.set_xticks([1,2])
    ax.set_xticklabels([f'Healthy\n(n={len(g_h)})', f'NSCLC\n(n={len(g_c)})'], fontsize=11)
    ax.set_title(f'{feat}\nMann–Whitney U, p = {p:.2e}  {pstar(p)}', fontweight='bold', fontsize=11)
    ax.set_ylabel(feat, fontsize=10); ax.grid(axis='y', alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')
    ax.spines[['top','right']].set_visible(False)

add_panel_labels(axes)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/eda_ch3_cancer_vs_healthy.png', dpi=150, bbox_inches='tight')
plt.show()

# ROC AUC for TF + คำนวณ OPT_TF, OPT_SENS, OPT_SPEC
y_roc  = (tp1_ichor['label']=='NSCLC').astype(int)
tf_roc = tp1_ichor['tumor_fraction'].fillna(0)
auc_tf = roc_auc_score(y_roc, tf_roc)
fpr_r, tpr_r, thr = roc_curve(y_roc, tf_roc)
opt_i   = np.argmax(tpr_r - fpr_r)
OPT_TF   = thr[opt_i]
OPT_SENS = tpr_r[opt_i]
OPT_SPEC = 1 - fpr_r[opt_i]
all_ichor['TF_status'] = (all_ichor['tumor_fraction'] > OPT_TF).astype(int)

print(f'tumor_fraction AUC (Healthy vs NSCLC): {auc_tf:.3f}')
print(f'Optimal threshold: {OPT_TF:.4f}  Sens={OPT_SENS:.3f}  Spec={OPT_SPEC:.3f}')

# %% [cell 11]
# ---
# ## Chapter 4 · Longitudinal Profiles by Response Group
#
# **Objective:** Compare ichorCNA features across PR/SD/PD at each timepoint.
#
# - 3 rows (TP1, TP2, TP3) × 3 columns (tumor_fraction, ploidy, gc_mad)
# - Kruskal–Wallis test for group differences per panel
# - Gray band = CRA healthy reference

# %% [cell 12]
can_ichor = (
    CRA_ichor[CRA_ichor['disease_status']=='lung cancer']
    .merge(outcome[['PID','CT']], left_on='pid', right_on='PID', how='left')
)
can_frag = (
    CRA_frag[CRA_frag['disease_status']=='lung cancer']
    .merge(outcome[['PID','CT']], left_on='pid', right_on='PID', how='left')
)
can_ichor['TF_status'] = (can_ichor['tumor_fraction'] > OPT_TF).astype(int)
hlt_tp1_ichor = CRA_ichor[CRA_ichor['disease_status']=='healthy']

fig, axes = plt.subplots(
    3, 3, figsize=(15, 14),
    gridspec_kw={'hspace': 0.52, 'wspace': 0.38}
)
fig.suptitle(
    'Longitudinal cfDNA Profiles by Treatment Response: Timepoint1, Timepoint2, and Timepoint3',
    fontsize=12, fontweight='bold', y=0.965
)

BOX_W  = 0.32
base_x = np.array([1, 2, 3])

for row, tp in enumerate([1, 2, 3]):
    tp_can = can_ichor[can_ichor['timepoint'] == tp]

    for col, feat in enumerate(ICHOR_FEATS):
        ax = axes[row, col]

        g_h      = hlt_tp1_ichor[feat].dropna().values
        grp_vals = [tp_can[tp_can['CT']==g][feat].dropna().values for g in GROUP_ORDER]

        # ── reference band ───────────────────────────────────────
        draw_ref_band(ax, g_h, color='#888888', alpha=0.18,
                      label=f'Healthy IQR (n={len(g_h)})')

        # ── boxplot ──────────────────────────────────────────────
        bp = ax.boxplot(grp_vals, positions=base_x, widths=BOX_W,
                        patch_artist=True, manage_ticks=False,
                        medianprops=dict(color='black', lw=2.2),
                        flierprops=dict(marker='o', ms=4, alpha=0.4))
        for box, g in zip(bp['boxes'], GROUP_ORDER):
            box.set_facecolor(PAL[g]); box.set_alpha(0.55)

        # ── jitter ───────────────────────────────────────────────
        np.random.seed(42)
        for xi, (vals, g) in enumerate(zip(grp_vals, GROUP_ORDER), 1):
            if len(vals) == 0: continue
            jx = xi + (np.random.rand(len(vals)) - 0.5) * BOX_W * 0.8
            ax.scatter(jx, vals, s=18, alpha=0.55, color=PAL[g],
                       edgecolors='none', zorder=4)

        all_data = np.concatenate(grp_vals + [g_h])
        ymin_d   = all_data.min()
        ymax_d   = all_data.max()
        yrange   = ymax_d - ymin_d if ymax_d != ymin_d else ymax_d * 0.1
        ax.set_ylim(ymin_d - yrange * 0.05,
                    ymax_d + yrange * 0.28)

        valid_v = [v for v in grp_vals if len(v) > 1]
        if len(valid_v) >= 2:
            _, p_kw = kruskal(*valid_v)
            kw_color = '#c0392b' if p_kw < 0.05 else '#666666'
            ax.text(0.5, 0.965,
                    f'KW p={p_kw:.3f}  {pstar(p_kw)}',
                    transform=ax.transAxes,
                    ha='center', va='top',
                    fontsize=8.5, fontweight='bold', color=kw_color,
                    bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                              alpha=0.80, edgecolor='#dddddd', linewidth=0.8))

        # ── x-axis labels ────────────────────────────────────────
        ax.set_xticks(base_x)
        ax.set_xticklabels(
            [f'{g}\n(n={len(tp_can[tp_can["CT"]==g])})' for g in GROUP_ORDER],
            fontsize=9.5)

        ax.set_ylabel(feat, fontsize=9.5)
        ax.grid(axis='y', alpha=0.3)
        ax.spines[['top','right']].set_visible(False)

        # feature title: top row only
        if row == 0:
            ax.set_title(feat, fontweight='bold', fontsize=11, pad=26)

        # legend: top-left panel only
        if row == 0 and col == 0:
            ax.legend(fontsize=7.5, loc='upper left', bbox_to_anchor=(0, 0.88))

    # ── Timepoint row label (right margin) ───────────────────────
    axes[row, 2].annotate(
        f'Timepoint{tp}',
        xy=(1.03, 0.5), xycoords='axes fraction',
        ha='left', va='center', fontsize=9, fontweight='bold',
        color='#2c3e50',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAF4FB',
                  edgecolor='#AED6F1', linewidth=0.8))

add_panel_labels(axes)
plt.savefig(f'{OUT}/eda_ch4_baseline_by_response.png', dpi=150, bbox_inches='tight')
plt.show()

print('Key observations:')
print('  • TP1 baseline: ทั้ง 3 features แทบไม่แยกกลุ่มได้ (KW ns) → ต้องใช้ trajectory')
print('  • tumor_fraction: PR ลดลงชัดที่ TP2/TP3; PD ยังคงสูง → discriminative ที่สุด')
print('  • ploidy & gc_mad: ไม่แยกกลุ่มได้ใน 3 timepoints')

# %% [cell 13]
# ---
# ## Chapter 5 · Longitudinal Trajectory
#
# **Objective:** Visualize individual patient trajectories and group-level trends over time.
#
# - **Left:** Spaghetti plot with group mean lines
# - **Right:** Boxplot distribution per timepoint with Kruskal–Wallis p-values
# - Features: `tumor_fraction`, `ploidy`

# %% [cell 14]
TPS      = [1, 2, 3]
TP_X     = {1: 1, 2: 3, 3: 5}
OFF_B    = {'PR': -0.38, 'SD': 0.0, 'PD': 0.38}
BW       = 0.32
hlt_ref  = {feat: CRA_ichor[CRA_ichor['disease_status']=='healthy'][feat].dropna().values
            for feat in ['tumor_fraction', 'ploidy']}

for feat, feat_label in [('tumor_fraction', 'Tumor Fraction'), ('ploidy', 'Ploidy')]:

    hlt_vals = hlt_ref[feat]
    all_vals = can_ichor[feat].dropna().values
    y_lo     = max(0, all_vals.min() - all_vals.max() * 0.03)
    y_hi     = all_vals.max() * 1.25

    fig = plt.figure(figsize=(16, 11))
    fig.suptitle('Longitudinal Trajectory by Treatment Response Group',
                 fontsize=12, fontweight='bold', y=1.00)

    gs = gridspec.GridSpec(2, 3, figure=fig,
                           hspace=0.50, wspace=0.30,
                           top=0.93, bottom=0.07,
                           left=0.07, right=0.97,
                           height_ratios=[1, 1])

    # ── ROW 0: spaghetti — แยกแต่ละ CT group ────────────────────────────────
    for col, ct_g in enumerate(GROUP_ORDER):
        ax = fig.add_subplot(gs[0, col])

        draw_ref_band(ax, hlt_vals, color='#888888', alpha=0.18,
                      label=f'Healthy IQR (n={len(hlt_vals)})')

        sub  = can_ichor[can_ichor['CT'] == ct_g].sort_values('timepoint')
        pids = sub['pid'].unique()

        for pid in pids:
            pt = sub[sub['pid'] == pid].sort_values('timepoint')
            if len(pt) < 2: continue
            ax.plot(pt['timepoint'].astype(int), pt[feat],
                    color=PAL[ct_g], alpha=0.50, lw=1.8,
                    marker='o', markersize=5,
                    markerfacecolor=PAL[ct_g],
                    markeredgecolor='white', markeredgewidth=0.6,
                    zorder=3)
            last = pt.iloc[-1]
            ax.annotate(pid.replace('2LB-', ''),
                        xy=(int(last['timepoint']), last[feat]),
                        xytext=(4, 1), textcoords='offset points',
                        fontsize=6, color=PAL[ct_g], alpha=0.65, va='bottom')

        ax.set_xticks(TPS)
        ax.set_xticklabels([TP_LABEL[t] for t in TPS], fontsize=9)
        ax.set_xlabel('Timepoint', fontsize=9)
        ax.set_ylabel(feat_label, fontsize=9)
        ax.set_ylim(y_lo, y_hi)
        ax.set_title(f'{ct_g}  (n = {len(pids)})',
                     fontweight='bold', fontsize=11, color=PAL[ct_g], pad=10)
        ax.grid(alpha=0.20)
        ax.spines[['top', 'right']].set_visible(False)
        if col == 0:
            ax.legend(fontsize=7.5, loc='upper right')

    # ── ROW 1: boxplot per TP (span all 3 cols) ──────────────────────────────
    ax_box = fig.add_subplot(gs[1, :])

    draw_ref_band(ax_box, hlt_vals, color='#888888', alpha=0.18,
                  label=f'Healthy IQR (n={len(hlt_vals)})')

    KW_p = {}
    np.random.seed(42)

    for tp in TPS:
        grp_vals_tp, pos_tp = [], []
        for g in GROUP_ORDER:
            v = can_ichor[(can_ichor['CT'] == g) &
                          (can_ichor['timepoint'] == tp)][feat].dropna().values
            grp_vals_tp.append(v)
            pos_tp.append(TP_X[tp] + OFF_B[g])

        bp = ax_box.boxplot(grp_vals_tp, positions=pos_tp, widths=BW,
                            patch_artist=True, manage_ticks=False, showfliers=False,
                            medianprops=dict(color='black', lw=2),
                            whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2))
        for box, g in zip(bp['boxes'], GROUP_ORDER):
            box.set_facecolor(PAL[g]); box.set_alpha(0.50)

        for xi, (v, g) in zip(pos_tp, zip(grp_vals_tp, GROUP_ORDER)):
            if len(v) == 0: continue
            jx = xi + (np.random.rand(len(v)) - 0.5) * BW * 0.7
            ax_box.scatter(jx, v, s=18, alpha=0.55, color=PAL[g],
                           edgecolors='none', zorder=4)

        valid_v = [v for v in grp_vals_tp if len(v) > 1]
        if len(valid_v) >= 2:
            _, p_kw = kruskal(*valid_v)
            KW_p[tp] = p_kw

    ax_box.set_ylim(y_lo, y_hi * 1.15)

    x_data_min = TP_X[1] + OFF_B['PR'] - BW
    x_data_max = TP_X[3] + OFF_B['PD'] + BW
    for tp in TPS:
        if tp not in KW_p: continue
        p_kw   = KW_p[tp]
        kw_col = '#c0392b' if p_kw < 0.05 else '#888888'
        x_ax   = (TP_X[tp] - x_data_min) / (x_data_max - x_data_min + 1e-9)
        x_ax   = max(0.05, min(0.95, x_ax))
        ax_box.text(x_ax, 0.965,
                    f'KW p={p_kw:.4f}  {pstar(p_kw)}',
                    transform=ax_box.transAxes,
                    ha='center', va='top',
                    fontsize=9, fontweight='bold', color=kw_col,
                    bbox=dict(boxstyle='round,pad=0.22', facecolor='white',
                              alpha=0.85, edgecolor='#dddddd', linewidth=0.8))

    ax_box.set_xticks([TP_X[t] for t in TPS])
    ax_box.set_xticklabels([TP_LABEL[t].replace('\n', ' ') for t in TPS], fontsize=10)
    ax_box.set_xlabel('Timepoint', fontsize=10)
    ax_box.set_ylabel(feat_label, fontsize=10)
    ax_box.set_title('Distribution by Timepoint  (* Kruskal–Wallis p < 0.05)',
                     fontweight='bold', fontsize=11)

    legend_h = [Patch(facecolor=PAL[g], alpha=0.65, label=g) for g in GROUP_ORDER]
    legend_h.append(Line2D([0], [0], color='#888888', lw=1.5, ls='--',
                           label=f'Healthy median (n={len(hlt_vals)})'))
    ax_box.legend(handles=legend_h, fontsize=9, loc='upper right')
    ax_box.grid(axis='y', alpha=0.25)
    ax_box.spines[['top', 'right']].set_visible(False)

    plt.tight_layout(pad=1.8)
    plt.savefig(f'{OUT}/eda_ch5_trajectory_{feat}.png', dpi=150, bbox_inches='tight')
    plt.show()

    print(f'\n── {feat_label} Kruskal-Wallis (PR vs SD vs PD) ──')

    tp_print_label = {
        1: 'Timepoint 1',
        2: 'Timepoint 2',
        3: 'Timepoint 3'
    }

    for tp in TPS:
        p = KW_p.get(tp, float('nan'))
        sig = '*** significant' if (not np.isnan(p) and p < 0.05) else '(ns)'
        print(f'  {tp_print_label[tp]:25s}: p = {p:.4f}  {sig}')

# %% [cell 15]
# ---
# ## Chapter 6 · Fragmentomics & End-Motif Analysis
#
# **Objective:** Examine fragment length features and end-motif patterns.
#
# - **6A:** 10 fragment length features — NSCLC vs Healthy (TP1)
# - **6B:** S/L Ratio longitudinal trajectory by response group
# - **6C:** End-motif PCA and top 20 differentially expressed motifs

# %% [cell 16]
FRAG_FEATS_ALL = ['frac_100_150','frac_151_180','frac_181_220','frac_221_300',
                  'short_to_long_100_150_over_151_220','skewness',
                  'kurtosis_excess','mean','median','short_fraction_s150']
FRAG_LABELS = {
    'frac_100_150':                       'Frac 100-150bp',
    'frac_151_180':                       'Frac 151-180bp',
    'frac_181_220':                       'Frac 181-220bp',
    'frac_221_300':                       'Frac 221-300bp',
    'short_to_long_100_150_over_151_220': 'S/L Ratio',
    'skewness':                           'Skewness',
    'kurtosis_excess':                    'Kurtosis',
    'mean':                               'Mean Length',
    'median':                             'Median Length',
    'short_fraction_s150':                'Short Frac (<150bp)',
}

# ── Healthy reference: TP1 (CRA+EGA combined, n=192) ─────────────────────────
hlt_frag_tp1 = all_frag[(all_frag['label']=='Healthy') & (all_frag['timepoint']==1)]

# ── CRA cancer with CT outcome ────────────────────────────────────────────────
can_frag_ct = (
    CRA_frag[CRA_frag['disease_status']=='lung cancer']
    .merge(outcome[['PID','CT']], left_on='pid', right_on='PID', how='left')
)

# ── x-position layout ─────────────────────────────────────────────────────────
# Healthy  |  TP1: PR SD PD  |  TP2: PR SD PD  |  TP3: PR SD PD
BW     = 0.22
GAP_TP = 0.55

def make_positions(bw=BW, gap=GAP_TP):
    pos  = {'Healthy': 0.0}
    base = 1.2
    for tp in [1, 2, 3]:
        for gi, g in enumerate(['PR','SD','PD']):
            pos[(tp, g)] = base + gi * (bw + 0.06)
        base += 3 * (bw + 0.06) + gap
    return pos

XPOS       = make_positions()
TP_CENTERS = {tp: np.mean([XPOS[(tp, g)] for g in GROUP_ORDER]) for tp in [1,2,3]}
HLT_X      = XPOS['Healthy']

# ── Figure: 10 rows × 1 col ───────────────────────────────────────────────────
n_feat = len(FRAG_FEATS_ALL)
fig, axes = plt.subplots(n_feat, 1, figsize=(14, n_feat * 3.2),
                          gridspec_kw={'hspace': 0.58})
fig.suptitle('Fragment Length Features across Timepoints: Healthy Controls and Response Groups',
             fontsize=12, fontweight='bold', y=0.92)

np.random.seed(42)

for row, feat in enumerate(FRAG_FEATS_ALL):
    ax       = axes[row]
    feat_lbl = FRAG_LABELS[feat]
    hlt_vals = hlt_frag_tp1[feat].dropna().values

    # shared y range
    all_v = can_frag_ct[feat].dropna().values
    combo = np.concatenate([all_v, hlt_vals])
    span  = combo.max() - combo.min() if combo.max() != combo.min() else abs(combo.max()) * 0.2
    ymin  = combo.min() - span * 0.05
    ymax  = combo.max() + span * 0.28

    # ── Healthy band + boxplot ────────────────────────────────────
    q1_h, q3_h = np.percentile(hlt_vals, [25, 75])
    ax.axhspan(q1_h, q3_h, color='#888888', alpha=0.13, zorder=0)
    ax.axhline(np.median(hlt_vals), color='#888888', lw=1.1, ls='--', alpha=0.65, zorder=1)

    bp_h = ax.boxplot([hlt_vals], positions=[HLT_X], widths=BW*1.1,
                      patch_artist=True, manage_ticks=False, showfliers=False,
                      medianprops=dict(color='black', lw=2.0),
                      whiskerprops=dict(lw=1.1), capprops=dict(lw=1.1))
    bp_h['boxes'][0].set_facecolor('#888888'); bp_h['boxes'][0].set_alpha(0.40)
    jx_h = HLT_X + (np.random.rand(len(hlt_vals)) - 0.5) * BW * 0.9
    ax.scatter(jx_h, hlt_vals, s=6, alpha=0.25, color='#888888', edgecolors='none', zorder=3)

    # ── Cancer boxes per TP ───────────────────────────────────────
    kw_per_tp = {}
    for tp in [1, 2, 3]:
        tp_data  = can_frag_ct[can_frag_ct['timepoint'] == tp]
        grp_vals = [tp_data[tp_data['CT']==g][feat].dropna().values for g in GROUP_ORDER]
        positions= [XPOS[(tp, g)] for g in GROUP_ORDER]

        bp = ax.boxplot(grp_vals, positions=positions, widths=BW,
                        patch_artist=True, manage_ticks=False, showfliers=False,
                        medianprops=dict(color='black', lw=2.0),
                        whiskerprops=dict(lw=1.1), capprops=dict(lw=1.1))
        for box, g in zip(bp['boxes'], GROUP_ORDER):
            box.set_facecolor(PAL[g]); box.set_alpha(0.58)

        for xi, (v, g) in zip(positions, zip(grp_vals, GROUP_ORDER)):
            if len(v) == 0: continue
            jx = xi + (np.random.rand(len(v)) - 0.5) * BW * 0.7
            ax.scatter(jx, v, s=14, alpha=0.55, color=PAL[g], edgecolors='none', zorder=4)

        valid_v = [v for v in grp_vals if len(v) > 1]
        if len(valid_v) >= 2:
            _, p_kw = kruskal(*valid_v)
            kw_per_tp[tp] = p_kw

    # ── cluster separator lines ───────────────────────────────────
    for sx in [0.75,
               XPOS[(1,'PD')] + BW/2 + 0.15,
               XPOS[(2,'PD')] + BW/2 + 0.15]:
        ax.axvline(sx, color='#cccccc', lw=0.9, ls=':', zorder=0)

    # ── KW labels: transAxes → ระนาบเดียวกันทุก panel ─────────────
    all_x     = list(XPOS.values())
    x_min_d   = min(all_x) - BW
    x_max_d   = max(all_x) + BW
    x_range_d = x_max_d - x_min_d

    for tp, p_kw in kw_per_tp.items():
        kw_col = '#c0392b' if p_kw < 0.05 else '#888888'
        x_ax   = (TP_CENTERS[tp] - x_min_d) / x_range_d
        ax.text(x_ax, 0.968,
                f'KW p={p_kw:.3f} {pstar(p_kw)}',
                transform=ax.transAxes, ha='center', va='top',
                fontsize=8, fontweight='bold', color=kw_col,
                bbox=dict(boxstyle='round,pad=0.18', facecolor='white',
                          alpha=0.88, edgecolor='#dddddd', linewidth=0.6))

    # ── x-axis labels ─────────────────────────────────────────────
    n_h  = len(hlt_vals)
    n_tp = {tp: {g: len(can_frag_ct[(can_frag_ct['CT']==g)&(can_frag_ct['timepoint']==tp)])
                 for g in GROUP_ORDER} for tp in [1,2,3]}
    ax.set_xticks([HLT_X] + [TP_CENTERS[tp] for tp in [1,2,3]])
    ax.set_xticklabels([
        f'Healthy\n(n={n_h})',
        f'Timepoint1',
        f'Timepoint2',
        f'Timepoint3',
    ], fontsize=9)

    # sub-label: PR / SD / PD color under each box
    for tp in [1, 2, 3]:
        for g in GROUP_ORDER:
            ax.text(XPOS[(tp, g)], ymin - span*0.03, g,
                    ha='center', va='top', fontsize=7.5,
                    color=PAL[g], fontweight='bold')

    ax.set_xlim(HLT_X - BW*2, XPOS[(3,'PD')] + BW*2)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(feat_lbl, fontsize=9.5)
    ax.set_title(feat_lbl, fontweight='bold', fontsize=10.5, pad=20, loc='left')
    ax.tick_params(axis='y', labelsize=8)
    ax.grid(axis='y', alpha=0.22, lw=0.6)
    ax.spines[['top','right']].set_visible(False)

# ── global legend ─────────────────────────────────────────────────
legend_h  = [Patch(facecolor=PAL[g], alpha=0.65, label=g) for g in GROUP_ORDER]
legend_h += [
    Patch(facecolor='#888888', alpha=0.40, label='Healthy (n=192)'),
    Line2D([0],[0], color='#888888', lw=1.2, ls='--', label='Healthy median'),
    Patch(facecolor='#888888', alpha=0.15, label='Healthy IQR'),
]
fig.legend(
    handles=legend_h,
    loc='lower center',
    ncol=6,
    fontsize=8.8,
    bbox_to_anchor=(0.5, 0.08),
    framealpha=0.85,
    handlelength=1.4,
    handletextpad=0.5,
    columnspacing=1.2,
    borderpad=0.45
)
add_panel_labels(axes)
plt.savefig(f'{OUT}/eda_ch6a_fraglen_cancer_vs_healthy.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [cell 17]
# 6B: Fragment length trajectory
hlt_frag_tp1 = CRA_frag[CRA_frag['disease_status']=='healthy']
feat_show = 'short_to_long_100_150_over_151_220'
feat_lbl  = 'S/L Ratio'

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle(
    'Longitudinal Trajectory by Treatment Response Group',
    fontsize=12, fontweight='bold', y=0.965
)

# -------------------------
# Panel A
# -------------------------
ax = axes[0]
hlt_ref_sl = hlt_frag_tp1[feat_show].dropna().values
draw_ref_band(ax, hlt_ref_sl, color='#888888', label=f'Healthy IQR (n={len(hlt_ref_sl)})')

for ct_g in GROUP_ORDER:
    sub = can_frag[can_frag['CT']==ct_g].sort_values('timepoint')
    for pid, grp in sub.groupby('pid'):
        g = grp.sort_values('timepoint')
        if len(g) < 2:
            continue
        ax.plot(
            g['timepoint'], g[feat_show],
            color=PAL[ct_g], alpha=0.28, lw=1.5,
            marker='o', markersize=3
        )

for ct_g in GROUP_ORDER:
    means = [
        can_frag[
            (can_frag['CT']==ct_g) & (can_frag['timepoint']==tp)
        ][feat_show].mean()
        for tp in TPS
    ]
    ax.plot(
        TPS, means,
        color=PAL[ct_g], lw=3, marker='o', ms=8,
        label=f'{ct_g} mean', zorder=5
    )

ax.set_xticks(TPS)
ax.set_xticklabels(['Timepoint1', 'Timepoint2', 'Timepoint3'], fontsize=10)
ax.set_ylabel(feat_lbl)
ax.set_title('Individual Trajectories with Group Mean', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(alpha=0.25)
ax.spines[['top','right']].set_visible(False)

# -------------------------
# Panel B
# -------------------------
ax = axes[1]
draw_ref_band(ax, hlt_ref_sl, color='#888888', label=f'Healthy IQR (n={len(hlt_ref_sl)})')

TP_X2 = {1:1, 2:3, 3:5}
OFF_6 = {'PR':-0.38, 'SD':0.0, 'PD':0.38}
BW6 = 0.32

np.random.seed(42)

for tp in TPS:
    gv = [
        can_frag[
            (can_frag['CT']==g) & (can_frag['timepoint']==tp)
        ][feat_show].dropna().values
        for g in GROUP_ORDER
    ]
    pos_v = [TP_X2[tp] + OFF_6[g] for g in GROUP_ORDER]

    bp = ax.boxplot(
        gv, positions=pos_v, widths=BW6,
        patch_artist=True, manage_ticks=False, showfliers=False,
        medianprops=dict(color='black', lw=2)
    )

    for box, g in zip(bp['boxes'], GROUP_ORDER):
        box.set_facecolor(PAL[g])
        box.set_alpha(0.50)

    for xi, (v, g) in zip(pos_v, zip(gv, GROUP_ORDER)):
        if len(v) == 0:
            continue
        jx = xi + (np.random.rand(len(v)) - 0.5) * BW6 * 0.7
        ax.scatter(jx, v, s=14, alpha=0.50, color=PAL[g],
                   edgecolors='none', zorder=4)

    valid = [v for v in gv if len(v) > 1]
    if len(valid) >= 2:
        _, p_kw = kruskal(*valid)
        ymax_tp = max(np.max(v) for v in gv if len(v) > 0)
        ax.text(
            TP_X2[tp], ymax_tp * 1.04, pstar(p_kw),
            ha='center', fontsize=9, fontweight='bold',
            color='#c0392b' if p_kw < 0.05 else '#888888'
        )

ax.set_xticks([TP_X2[t] for t in TPS])
ax.set_xticklabels(['Timepoint1', 'Timepoint2', 'Timepoint3'], fontsize=10)
ax.set_ylabel(feat_lbl)
ax.set_title('Boxplot per Timepoint', fontweight='bold')

legend_h = [Patch(facecolor=PAL[g], alpha=0.6, label=g) for g in GROUP_ORDER]
legend_h.append(Line2D([0], [0], color='#888888', lw=1.5, ls='--', label='Healthy median'))

ax.legend(handles=legend_h, fontsize=9)
ax.grid(axis='y', alpha=0.25)
ax.spines[['top','right']].set_visible(False)

add_panel_labels(axes)
plt.tight_layout(rect=[0.02, 0.03, 0.98, 0.93], pad=1.8)
plt.savefig(f'{OUT}/eda_ch6b_fraglen_trajectory.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [cell 18]
# 6C: End motif PCA + Top differential motifs
motif_tp1 = all_motif[all_motif['timepoint']==1].dropna(subset=MOTIF_COLS)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('Cell-Free DNA End-Motif Profiles: NSCLC versus Healthy Controls',
             fontsize=13, fontweight='bold', y=1.01)

ax = axes[0]
Xm = StandardScaler().fit_transform(motif_tp1[MOTIF_COLS].values)
pca2 = PCA(n_components=2); pcs = pca2.fit_transform(Xm)
for lbl, col in [('Healthy',PAL['Healthy']),('NSCLC',PAL['NSCLC'])]:
    mask = motif_tp1['label'].values == lbl
    ax.scatter(pcs[mask,0], pcs[mask,1], c=col, label=lbl,
               alpha=0.7, s=45, edgecolors='white', lw=0.5)
ax.set_xlabel(f'PC1 ({pca2.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca2.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('Principal Component Analysis of End-Motif Frequencies', fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.25); ax.spines[['top','right']].set_visible(False)

ax = axes[1]
nm_ = motif_tp1[motif_tp1['label']=='NSCLC'][MOTIF_COLS]
hm_ = motif_tp1[motif_tp1['label']=='Healthy'][MOTIF_COLS]
mpv = []
for m in MOTIF_COLS:
    _, p = mannwhitneyu(nm_[m].dropna(), hm_[m].dropna(), alternative='two-sided')
    mpv.append((m, p, nm_[m].median()-hm_[m].median()))
mpv_df = pd.DataFrame(mpv, columns=['motif','p','diff']).sort_values('p').head(20)
cols_b = [PAL['NSCLC'] if d > 0 else PAL['Healthy'] for d in mpv_df['diff']]
ax.barh(mpv_df['motif'][::-1], mpv_df['diff'][::-1],
        color=cols_b[::-1], alpha=0.80, edgecolor='white')
ax.axvline(0, color='black', lw=1)
ax.set_title('Top 20 Differentially Expressed End Motifs\n(NSCLC − Healthy Controls)', fontweight='bold')
ax.set_xlabel('Median Difference'); ax.grid(axis='x', alpha=0.3)
ax.spines[['top','right']].set_visible(False)
add_panel_labels(axes)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/eda_ch6c_endmotif.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [cell 19]
# ---
# ## Chapter 7 · Feature Discriminability Summary
#
# **Objective:** Identify which features best separate PR/SD/PD across timepoints.
#
# - Heatmap: Kruskal–Wallis p-value for each feature × timepoint combination
# - 🔴 Red = low p (significant), 🟢 Green = high p (not significant)
# - Features marked `*` are significant at p < 0.05

# %% [cell 20]
FEATS_ALL = {
    'ichorCNA': ['tumor_fraction','ploidy','gc_mad'],
    'fraglen':  ['short_to_long_100_150_over_151_220','short_fraction_s150',
                 'mean','skewness','frac_151_180','frac_181_220'],
}

summary_rows = []
for category, feats in FEATS_ALL.items():
    src = can_ichor if category == 'ichorCNA' else can_frag
    for feat in feats:
        for tp in [1, 2, 3]:
            vals = [
                src[(src['CT'] == g) & (src['timepoint'] == tp)][feat].dropna().values
                for g in GROUP_ORDER
            ]
            p = kruskal(*vals)[1] if all(len(v) > 1 for v in vals) else float('nan')
            summary_rows.append({
                'category': category,
                'feature': feat,
                'TP': tp,
                'p_kw': p
            })

kw_summary = pd.DataFrame(summary_rows)
kw_pivot = kw_summary.pivot_table(
    index=['category', 'feature'],
    columns='TP',
    values='p_kw'
)
kw_pivot.columns = ['TP1_p', 'TP2_p', 'TP3_p']

# ==============================
# Heatmap
# ==============================
fig, ax = plt.subplots(figsize=(9.2, 6.2))

fig.suptitle(
    'Feature Discriminability across Timepoints (Kruskal-Wallis p-value)',
    fontsize=13,
    fontweight='bold',
    y=1.01
)

heat_data = kw_pivot.reset_index()
feat_labels = [
    f"{r['category']}\n{r['feature']}"
    for _, r in heat_data.iterrows()
]

mat = heat_data[['TP1_p', 'TP2_p', 'TP3_p']].values

im = ax.imshow(
    mat,
    cmap=plt.cm.RdYlGn_r,
    vmin=0,
    vmax=0.15,
    aspect='auto'
)

ax.set_xticks([0, 1, 2])

ax.set_xticklabels(
    [TP_LABEL[t].replace('\n', ' ') for t in [1, 2, 3]],
    fontsize=9.5,
    rotation=15,
    ha='right',
    rotation_mode='anchor'
)

ax.tick_params(axis='x', pad=8)

ax.set_yticks(range(len(feat_labels)))
ax.set_yticklabels(feat_labels, fontsize=9)

for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat[i, j]
        text = f'{v:.3f}' if not np.isnan(v) else 'NA'
        star = ' *' if (not np.isnan(v) and v < 0.05) else ''

        ax.text(
            j, i, text + star,
            ha='center',
            va='center',
            fontsize=8.5,
            fontweight='bold' if (not np.isnan(v) and v < 0.05) else 'normal',
            color='white' if (not np.isnan(v) and v < 0.03) else 'black'
        )

plt.colorbar(
    im,
    ax=ax,
    shrink=0.7,
    label='KW p-value'
)

ax.axhline(2.5, color='white', lw=2)

add_panel_labels(ax)
plt.tight_layout(pad=1.8)
plt.subplots_adjust(bottom=0.16)

plt.savefig(
    f'{OUT}/eda_ch7_feature_summary_heatmap.png',
    dpi=150,
    bbox_inches='tight'
)

plt.show()

display(
    kw_pivot.style
    .format('{:.4f}')
    .background_gradient(cmap='RdYlGn_r', vmin=0, vmax=0.15)
    .set_caption('Kruskal-Wallis p-value: PR vs SD vs PD')
)

# %% [cell 21]
# ── Display Helper: Styled Summary Table ──────────────────────────────────
def styled_summary(df, title='Summary', gradient_cols=None, vmin=0, vmax=1):
    """Display a styled summary DataFrame with optional gradient coloring."""
    styler = df.style.set_caption(title)
    if gradient_cols:
        styler = styler.background_gradient(
            subset=gradient_cols, cmap='RdYlGn', vmin=vmin, vmax=vmax
        )
    return styler

def print_section(title, char='═', width=65):
    """Print a formatted section header."""
    print(f'\n{char * width}')
    print(f'  {title}')
    print(f'{char * width}')

print('✓ Display helpers loaded')

# %% [cell 22]
# ---
# # Part 2 · Feature Engineering
#
# Build a **patient-level feature matrix** (1 row per patient) combining all data modalities.
#
# | Feature Group | Examples | Count |
# |--------------|---------|-------|
# | Raw per timepoint | `tumor_fraction_t1`, `skewness_t2` | 3 × TP × features |
# | Absolute delta | `d_tumor_fraction_dt2` | Δ T0→W3, T0→W12, W3→W12 |
# | % change | `pct_tumor_fraction_dt2` | Same pairs |
# | TF binary flags | `TF_pos_t1`, `TF_clear_t2`, `TF_rebound_t3` | 5 flags |
# | End-motif PCA | `mpc1_t1 ... mpc10_t3` | 10 PCs × 3 TP |
# | Cancer probability (LOOCV) | `cp_t1/t2/t3`, `dcp_dt2/dt3` | 5 features |
#
# > ⚠️ **Anti-leakage:** End-motif PCA is fitted on Healthy samples only.

# %% [cell 23]
from sklearn.impute import SimpleImputer

# Healthy-only PCA (anti-leakage)
hlt_motif = all_motif[all_motif['label']=='Healthy']
Xm_hlt    = hlt_motif[MOTIF_COLS].fillna(0).values
sc_m      = StandardScaler().fit(Xm_hlt)
pca10     = PCA(n_components=10, random_state=SEED).fit(sc_m.transform(Xm_hlt))
print(f'PCA on {len(Xm_hlt)} Healthy  PC1={pca10.explained_variance_ratio_[0]*100:.1f}%')

def motif_pca_row(pid, tp):
    r = CRA_motif[(CRA_motif['pid']==pid)&(CRA_motif['timepoint']==tp)]
    if len(r)==0: return {f'mpc{k}': np.nan for k in range(1,11)}
    pc = pca10.transform(sc_m.transform(r[MOTIF_COLS].fillna(0).values))
    return {f'mpc{k}': pc[0][k-1] for k in range(1,11)}

cancer_pids = [p for p in CRA_ichor[CRA_ichor['disease_status']=='lung cancer']['pid'].unique()
               if p in outcome['PID'].values]
print(f'Cancer patients with outcome: {len(cancer_pids)}')

DELTA_PAIRS = [(1,2,'dt2'),(1,3,'dt3'),(2,3,'dt23')]
rows = []

for pid in cancer_pids:
    ct  = outcome[outcome['PID']==pid]['CT'].values[0]
    row = {'pid':pid,'CT':ct}

    for tp in [1,2,3]:
        ri = CRA_ichor[(CRA_ichor['pid']==pid)&(CRA_ichor['timepoint']==tp)]
        rf = CRA_frag[ (CRA_frag['pid']==pid) &(CRA_frag['timepoint']==tp)]
        for f in ICHOR_FEATS:
            row[f'{f}_t{tp}'] = ri[f].values[0] if len(ri)>0 else np.nan
        tfv = ri['tumor_fraction'].values[0] if len(ri)>0 else np.nan
        row[f'TF_pos_t{tp}'] = int(tfv>OPT_TF) if not np.isnan(tfv) else np.nan
        for f in FRAG_FEATS:
            row[f'{f}_t{tp}'] = rf[f].values[0] if len(rf)>0 else np.nan
        for k,v in motif_pca_row(pid,tp).items():
            row[f'{k}_t{tp}'] = v

    for f in ICHOR_FEATS+FRAG_FEATS:
        for t1,t2,tag in DELTA_PAIRS:
            v1=row.get(f'{f}_t{t1}',np.nan); v2=row.get(f'{f}_t{t2}',np.nan)
            ok = not(np.isnan(v1) or np.isnan(v2))
            row[f'd_{f}_{tag}']   = (v2-v1) if ok else np.nan
            row[f'pct_{f}_{tag}'] = ((v2-v1)/abs(v1)*100 if(ok and v1!=0) else np.nan)

    row['TF_clear_t2']   = int(row.get('TF_pos_t1',0)==1 and row.get('TF_pos_t2',1)==0)
    row['TF_clear_t3']   = int(row.get('TF_pos_t1',0)==1 and row.get('TF_pos_t3',1)==0)
    row['TF_rebound_t3'] = int(row.get('TF_pos_t2',0)==0 and row.get('TF_pos_t3',1)==1)
    rows.append(row)

feat_df = pd.DataFrame(rows)
all_fc  = [c for c in feat_df.columns if c not in ['pid','CT']]
print(f'Feature matrix: {feat_df.shape}  response={feat_df["CT"].value_counts().to_dict()}')

# %% [cell 24]
# ---
# # Part 3 · Stage 1: Cancer vs Healthy Classification
#
# **Task:** Binary classification — NSCLC vs Healthy  
# **Data:** CRA + EGA combined (all timepoints → TP1 features only)  
# **Method:** Stratified 5-Fold CV with SMOTE oversampling  
# **Models:** LR, RF, SVM, GB, LGBM
#
# | Metric | Description |
# |--------|-------------|
# | AUC | Area under ROC curve |
# | Sens | Sensitivity (Recall for Cancer) |
# | Spec | Specificity |
# | F2 | F-beta with β=2 (emphasizes recall) |
# | BalAcc | Balanced Accuracy |
#
# > 💡 Best model (by BalAcc) is used to generate `cp_t1/t2/t3` cancer probabilities via LOOCV.

# %% [cell 25]
# ── Seed Reset Checkpoint ─────────────────────────────────────────────────────
set_global_seed(SEED)
print(f'✓ All random seeds reset to {SEED} — results will be reproducible')

# %% [cell 26]
from sklearn.metrics import fbeta_score, balanced_accuracy_score

MPC_COLS = [f'mpc{k}' for k in range(1,11)]

def build_s1(ichor, frag, motif):
    base = (ichor[ichor['timepoint']==1][['sample_key','disease_status','pid']+ICHOR_FEATS]
            .merge(frag[frag['timepoint']==1][['sample_key']+FRAG_FEATS],on='sample_key',how='inner'))
    mt = motif[motif['timepoint']==1].copy()
    Xm = mt[MOTIF_COLS].fillna(0).values
    pcs = pca10.transform(sc_m.transform(Xm))
    mpc_df = pd.DataFrame(pcs, columns=MPC_COLS)
    mpc_df['sample_key'] = mt['sample_key'].values
    return base.merge(mpc_df, on='sample_key', how='left')

s1  = pd.concat([build_s1(CRA_ichor,CRA_frag,CRA_motif),
                 build_s1(EGA_ichor,EGA_frag,EGA_motif)], ignore_index=True)
S1F = ICHOR_FEATS+FRAG_FEATS+MPC_COLS
X1  = s1[S1F].values
y1  = (s1['disease_status']=='lung cancer').astype(int).values
print(f'S1: n={len(s1)} Cancer={y1.sum()} Healthy={(y1==0).sum()} ratio={(y1==0).sum()/y1.sum():.1f}:1')

def smote_oversample(X,y,k=5,rs=42):
    rng=np.random.RandomState(rs); mn=np.where(y==1)[0]
    n=int((y==0).sum())-int((y==1).sum()); syn=[]
    for _ in range(n):
        idx=rng.choice(mn); d=np.linalg.norm(X[mn]-X[idx],axis=1)
        nb_=mn[np.argsort(d)[1:k+1]]; lam=rng.uniform(0,1)
        syn.append(X[idx]+lam*(X[rng.choice(nb_)]-X[idx]))
    Xs=np.array(syn)
    return np.vstack([X,Xs]), np.hstack([y,np.ones(len(Xs),dtype=int)])

def bsw(y_arr):
    cls,cnt=np.unique(y_arr,return_counts=True); total=len(y_arr); w=np.ones(len(y_arr))
    for c,n in zip(cls,cnt): w[y_arr==c]=total/(len(cls)*n)
    return w

s1_clfs = {
    'LR':   LogisticRegression(C=1.0, max_iter=1000, random_state=42),
    'RF':   RandomForestClassifier(n_estimators=200, random_state=42),
    'SVM':  SVC(probability=True, random_state=42),
    'GB':   GradientBoostingClassifier(n_estimators=100, random_state=42),
    'LGBM': lgb.LGBMClassifier(
                n_estimators=200, learning_rate=0.05, max_depth=4,
                num_leaves=15, min_child_samples=3,
                class_weight='balanced', random_state=42,
                verbosity=-1, n_jobs=-1),
}

IMP_S1=SimpleImputer(strategy='median'); SC_S1=StandardScaler()
cv5=StratifiedKFold(5,shuffle=True,random_state=SEED); s1_res={}
set_global_seed(SEED)  # reset before CV loop

print(f"\n{'Model':5} {'AUC':>6} {'Sens':>6} {'Spec':>6} {'F2':>6} {'BalAcc':>7}  thr")
print('-'*50)
for nm,clf in s1_clfs.items():
    yp=np.zeros(len(y1))
    for tr,te in cv5.split(X1,y1):
        Xi=IMP_S1.fit_transform(X1[tr]); Xe=IMP_S1.transform(X1[te])
        Xs_=SC_S1.fit_transform(Xi);     Xe_=SC_S1.transform(Xe)
        Xsm,ysm=smote_oversample(Xs_,y1[tr],k=5,rs=42)
        if nm in ('GB','LGBM'): clf.fit(Xsm,ysm,sample_weight=bsw(ysm))
        else:                    clf.fit(Xsm,ysm)
        yp[te]=clf.predict_proba(Xe_)[:,1]
    fpr,tpr,thrs=roc_curve(y1,yp); thr=float(thrs[np.argmax(tpr-fpr)]); yd=(yp>thr).astype(int)
    s1_res[nm]={'AUC':roc_auc_score(y1,yp),'Sens':recall_score(y1,yd,zero_division=0),
                'Spec':recall_score(1-y1,1-yd,zero_division=0),'F2':fbeta_score(y1,yd,beta=2,zero_division=0),
                'BalAcc':balanced_accuracy_score(y1,yd),'thr':thr,'yp':yp,'yd':yd}
    r=s1_res[nm]
    print(f"{nm:5} {r['AUC']:6.3f} {r['Sens']:6.3f} {r['Spec']:6.3f} {r['F2']:6.3f} {r['BalAcc']:7.3f}  {thr:.3f}")

best_s1=max(s1_res,key=lambda k:s1_res[k]['BalAcc'])
print(f'\nBest (BalAcc): {best_s1}  AUC={s1_res[best_s1]["AUC"]:.3f}')
S1_THR=s1_res[best_s1]['thr']
X1_imp=IMP_S1.fit_transform(X1); X1_sc=SC_S1.fit_transform(X1_imp)
X1_sm,y1_sm=smote_oversample(X1_sc,y1,k=5,rs=42)
best_clf_s1=s1_clfs[best_s1]
if best_s1=='GB': best_clf_s1.fit(X1_sm,y1_sm,sample_weight=bsw(y1_sm))
else:             best_clf_s1.fit(X1_sm,y1_sm)

# %% [cell 27]
# ### Stage 1 · ROC + Confusion Matrices

# %% [cell 28]
colors_s1 = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']  # +1 สีสำหรับ LGBM
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
fig.patch.set_facecolor('white')
fig.suptitle('Stage 1 Classification: NSCLC versus Healthy Controls',
             fontsize=14, fontweight='bold', y=1.01)

ax = axes[0]
for (nm, r), col in zip(s1_res.items(), colors_s1):
    fp, tp_c, _ = roc_curve(y1, r['yp'])
    ax.plot(fp, tp_c, color=col, lw=2, label=f"{nm} {r['AUC']:.3f}")
    oi = np.argmax(tp_c - fp)
    ax.scatter(fp[oi], tp_c[oi], color=col, s=60, zorder=5, edgecolors='white', lw=1.2)
ax.plot([0,1],[0,1],'k--',lw=1)
ax.set_xlabel('1 − Specificity  (False Positive Rate)')
ax.set_ylabel('Sensitivity  (True Positive Rate)')
ax.set_title('ROC Curve  (dot = Youden Index threshold)', fontweight='bold')
ax.legend(fontsize=9, loc='lower right')
ax.grid(alpha=0.3); ax.spines[['top','right']].set_visible(False)

ax = axes[1]
rf_s1 = RandomForestClassifier(n_estimators=200, random_state=42).fit(X1_sm, y1_sm)
imp_s1 = pd.Series(rf_s1.feature_importances_, index=S1F).sort_values(ascending=False).head(12)
ax.barh(imp_s1.index[::-1], imp_s1.values[::-1], color='#888888', alpha=0.85, edgecolor='white')
ax.set_title('Top 12 Feature Importances\n(Random Forest + SMOTE)', fontweight='bold', fontsize=11)
ax.set_xlabel('Importance'); ax.grid(axis='x', alpha=0.3); ax.spines[['top','right']].set_visible(False)
add_panel_labels(axes)
plt.tight_layout()
plt.savefig(f'{OUT}/fig_s1_roc.png', dpi=150, bbox_inches='tight')
plt.show()


fig2, axes2 = plt.subplots(1, 5, figsize=(25, 4.5))
fig2.suptitle('Confusion Matrices: Stage 1 Classification at Youden Threshold',
              fontsize=13, fontweight='bold')
for ax, (nm, r), col in zip(axes2, s1_res.items(), colors_s1):
    ConfusionMatrixDisplay(confusion_matrix(y1, r['yd']),
                           display_labels=['Healthy', 'Cancer']).plot(
        ax=ax, colorbar=False, cmap='Blues')
    ax.set_title(f"{nm}\nAUC={r['AUC']:.3f} thr={r['thr']:.2f}\n"
                 f"Sens={r['Sens']:.2f} Spec={r['Spec']:.2f}",
                 fontweight='bold', fontsize=9)
plt.tight_layout()
plt.savefig(f'{OUT}/fig_s1_confusion.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [cell 29]
# ### Stage 1 · SHAP Feature Explanation (Best Model)

# %% [cell 30]
# ── Stage 1 · SHAP Explanation ───────────────────────────────────────────
# Train best Stage-1 model on full dataset for global SHAP explanation

import shap
import lightgbm as lgb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# -----------------------------
# Prepare data
# -----------------------------
imp_sh = SimpleImputer(strategy='median')
sc_sh  = StandardScaler()

X1_imp = sc_sh.fit_transform(imp_sh.fit_transform(X1))
X1_sm, y1_sm = smote_oversample(X1_imp, y1, k=5, rs=42)
sw1_sm = bsw(y1_sm)

def _make_s1_shap_clf():
    if best_s1 == 'LR':
        return LogisticRegression(C=1.0, max_iter=1000)
    if best_s1 == 'RF':
        return RandomForestClassifier(n_estimators=200, random_state=42)
    if best_s1 == 'SVM':
        return SVC(probability=True, random_state=SEED)
    if best_s1 == 'LGBM':
        return lgb.LGBMClassifier(
            n_estimators=200,
            learning_rate=0.05,
            max_depth=4,
            num_leaves=15,
            min_child_samples=3,
            class_weight='balanced',
            random_state=42,
            verbosity=-1
        )
    return GradientBoostingClassifier(n_estimators=100, random_state=42)

clf_sh = _make_s1_shap_clf()
if best_s1 in ('GB', 'LGBM'):
    clf_sh.fit(X1_sm, y1_sm, sample_weight=sw1_sm)
else:
    clf_sh.fit(X1_sm, y1_sm)

# -----------------------------
# Build SHAP explainer
# -----------------------------
if best_s1 == 'LR':
    explainer = shap.LinearExplainer(
        clf_sh, X1_imp, feature_perturbation='interventional'
    )
    shap_vals = explainer.shap_values(X1_imp)

elif best_s1 in ('RF', 'GB', 'LGBM'):
    explainer = shap.TreeExplainer(clf_sh)
    sv = explainer.shap_values(X1_imp)

    if isinstance(sv, list):
        shap_vals = sv[1]         
    elif hasattr(sv, 'ndim') and sv.ndim == 3:
        shap_vals = sv[:, :, 1]
    else:
        shap_vals = sv

else:  
    bg = shap.kmeans(X1_imp, 20)
    explainer = shap.KernelExplainer(clf_sh.predict_proba, bg)
    sv_svm = explainer.shap_values(X1_imp, nsamples=100)

    if isinstance(sv_svm, list):
        shap_vals = sv_svm[1]
    elif hasattr(sv_svm, 'ndim') and sv_svm.ndim == 3:
        shap_vals = sv_svm[:, :, 1]
    else:
        shap_vals = sv_svm

print(f'[SHAP S1] best_s1={best_s1}  shap_vals shape={shap_vals.shape}')

# -----------------------------
# Feature display names
# -----------------------------
DISPLAY_NAMES = {
    'tumor_fraction': 'tumor_fraction',
    'gc_mad': 'gc_mad',
    'ploidy': 'ploidy',
    'short_fraction_s150': 'short_frac_s150',
    'short_to_long_100_150_over_151_220': 'S/L ratio',
    'frac_100_150': 'frac_100_150',
    'frac_151_180': 'frac_151_180',
    'frac_181_220': 'frac_181_220',
    'frac_221_300': 'frac_221_300',
    'mean': 'mean_len',
    'median': 'median_len',
    'skewness': 'skewness',
    'kurtosis_excess': 'kurtosis',
}

# -----------------------------
# Top 15 features by mean |SHAP|
# -----------------------------
mean_abs_all = np.abs(shap_vals).mean(axis=0)
top_idx = np.argsort(mean_abs_all)[-15:][::-1]

top_features = [S1F[i] for i in top_idx]
top_labels   = [DISPLAY_NAMES.get(f, f) for f in top_features]
top_shap     = shap_vals[:, top_idx]
top_X        = X1_imp[:, top_idx]
top_mean_abs = mean_abs_all[top_idx]


shap_cmap = LinearSegmentedColormap.from_list(
    'shap_blue_pink',
    ['#1E88E5', '#FF0052']
)

# -----------------------------
# Create figure layout
# -----------------------------
fig = plt.figure(figsize=(14.5, 7.8))
gs = fig.add_gridspec(1, 3, width_ratios=[2.35, 0.10, 1.75], wspace=0.55)

ax_bee = fig.add_subplot(gs[0, 0])
cax    = fig.add_subplot(gs[0, 1])
ax_bar = fig.add_subplot(gs[0, 2])

fig.suptitle(
    f'SHAP Feature Explanation — Stage 1 Model ({best_s1}): NSCLC versus Healthy',
    fontsize=14, fontweight='bold', y=0.99
)

# -----------------------------
# Left panel: custom beeswarm
# -----------------------------
rng = np.random.default_rng(42)
n_top = len(top_idx)

for row in range(n_top):
    sv = top_shap[:, row]
    fv = top_X[:, row]

    lo = np.nanpercentile(fv, 5)
    hi = np.nanpercentile(fv, 95)
    if hi > lo:
        fv_norm = np.clip((fv - lo) / (hi - lo), 0, 1)
    else:
        fv_norm = np.full_like(fv, 0.5, dtype=float)

    y = row + rng.uniform(-0.22, 0.22, size=len(sv))

    ax_bee.scatter(
        sv, y,
        c=fv_norm, cmap=shap_cmap,
        s=16, alpha=0.95, linewidths=0, zorder=3
    )

ax_bee.axvline(0, color='#888888', lw=1.2, ls='--', zorder=1)
ax_bee.set_yticks(range(n_top))
ax_bee.set_yticklabels(top_labels, fontsize=10.5)
ax_bee.invert_yaxis()

ax_bee.set_title('SHAP Beeswarm Plot — Top 15 Features',
                 fontsize=13, fontweight='bold', pad=8)
ax_bee.set_xlabel('SHAP value (impact on model output)', fontsize=11)
ax_bee.grid(axis='y', alpha=0.18, linestyle=':', linewidth=0.8)
ax_bee.grid(axis='x', alpha=0.15, linestyle='--', linewidth=0.8)
ax_bee.spines[['top', 'right', 'left']].set_visible(False)
ax_bee.tick_params(axis='y', pad=6)
ax_bee.tick_params(axis='x', labelsize=10)

# -----------------------------
# Middle: colorbar
# -----------------------------
sm = plt.cm.ScalarMappable(cmap=shap_cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['Low', 'High'])
cbar.set_label('')
cbar.ax.set_title('Feature value', fontsize=11, pad=8)

cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=10)

# -----------------------------
# Right panel: mean |SHAP| bar plot
# -----------------------------
ypos = np.arange(n_top)

ax_bar.barh(
    ypos, top_mean_abs,
    color='#1E88E5',
    height=0.68
)

ax_bar.set_yticks(ypos)
ax_bar.set_yticklabels(top_labels, fontsize=10.5)
ax_bar.invert_yaxis()

ax_bar.set_title('Mean |SHAP| Feature Importance',
                 fontsize=13, fontweight='bold', pad=8)
ax_bar.set_xlabel('mean |SHAP value|', fontsize=11)
ax_bar.grid(axis='x', alpha=0.18, linestyle='--', linewidth=0.8)
ax_bar.spines[['top', 'right', 'left']].set_visible(False)
ax_bar.tick_params(axis='y', pad=4)
ax_bar.tick_params(axis='x', labelsize=10)

# -----------------------------
# Final layout
# -----------------------------
fig.subplots_adjust(left=0.24, right=0.97, top=0.92, bottom=0.11)
plt.savefig(f'{OUT}fig_shap_stage1.png', dpi=150, bbox_inches='tight')
plt.show()

print('✓ fig_shap_stage1.png')

# -----------------------------
# SHAP importance table
# -----------------------------
shap_imp_s1 = pd.Series(
    np.abs(shap_vals).mean(axis=0),
    index=S1F
).sort_values(ascending=False)

print('Top-10 SHAP features (Stage 1):')
print(shap_imp_s1.head(10).to_string())

# %% [cell 31]
# ### Stage 1 · LOOCV Cancer Probability (anti-leakage)

# %% [cell 32]
s1_cancer=s1[s1['disease_status']=='lung cancer'].copy()
s1_healthy=s1[s1['disease_status']=='healthy'].copy()

def make_s1_clf():
    if best_s1=='LR':  return LogisticRegression(C=1.0,max_iter=1000)
    if best_s1=='RF':  return RandomForestClassifier(n_estimators=200,random_state=42)
    if best_s1=='SVM': return SVC(probability=True, random_state=SEED)
    if best_s1=='LGBM': return lgb.LGBMClassifier(
        n_estimators=200,learning_rate=0.05,max_depth=4,num_leaves=15,
        min_child_samples=3,class_weight='balanced',random_state=42,
        verbosity=-1,n_jobs=-1)
    return GradientBoostingClassifier(n_estimators=100,random_state=42)

def get_cp_loo(pid,tp):
    ri=CRA_ichor[(CRA_ichor['pid']==pid)&(CRA_ichor['timepoint']==tp)&(CRA_ichor['disease_status']=='lung cancer')]
    rf=CRA_frag[ (CRA_frag['pid']==pid) &(CRA_frag['timepoint']==tp) &(CRA_frag['disease_status']=='lung cancer')]
    rm=CRA_motif[(CRA_motif['pid']==pid)&(CRA_motif['timepoint']==tp)]
    if len(ri)==0 or len(rf)==0: return np.nan
    rd={**{f:ri[f].values[0] for f in ICHOR_FEATS}, **{f:rf[f].values[0] for f in FRAG_FEATS}}
    if len(rm)>0:
        pcs=pca10.transform(sc_m.transform(rm[MOTIF_COLS].fillna(0).values))[0]
        rd.update({f'mpc{k+1}':pcs[k] for k in range(10)})
    else: rd.update({f'mpc{k}':0.0 for k in range(1,11)})
    xd=pd.DataFrame([rd])[S1F]
    s1_tr=pd.concat([s1_healthy,s1_cancer[s1_cancer['pid']!=pid]],ignore_index=True)
    Xtr=s1_tr[S1F].values; ytr=(s1_tr['disease_status']=='lung cancer').astype(int).values
    imp_=SimpleImputer(strategy='median'); sc_=StandardScaler()
    Xt=sc_.fit_transform(imp_.fit_transform(Xtr)); xds=sc_.transform(imp_.transform(xd.values))
    Xs,ys=smote_oversample(Xt,ytr,k=5,rs=42)
    clf_=make_s1_clf()
    if best_s1 in ('GB','LGBM'): clf_.fit(Xs,ys,sample_weight=bsw(ys))
    else:                         clf_.fit(Xs,ys)
    return float(clf_.predict_proba(xds)[0,1])

set_global_seed(SEED)  # reset before LOOCV
print(f'LOOCV cp: {len(cancer_pids)} patients x 3 timepoints ...')
for tp in [1,2,3]:
    feat_df[f'cp_t{tp}']=feat_df['pid'].apply(lambda p: get_cp_loo(p,tp))
    print(f'  tp={tp} done  median={feat_df["cp_t"+str(tp)].median():.3f}')
feat_df['dcp_dt2']=feat_df['cp_t2']-feat_df['cp_t1']
feat_df['dcp_dt3']=feat_df['cp_t3']-feat_df['cp_t1']
print(f'Total features: {feat_df.shape[1]-2}')
display(feat_df[['pid','CT','cp_t1','cp_t2','cp_t3','dcp_dt2','dcp_dt3']].head())

# %% [cell 33]
# ---
# # Part 4 · Stage 2: Treatment Response Prediction
#
# **Task:** PR vs SD+PD classification  
# **Method:** LOOCV (Leave-One-Out CV) across 5 scenarios × 4 feature types × 5 models
#
# **Scenarios:**
# | Scenario | Data Used | Clinical Question |
# |----------|-----------|------------------|
# | `baseline` | T0 only | Can we predict before treatment? |
# | `early` | W3-4 only | Early response signal |
# | `late` | W12 only | Late response signal |
# | `change` | Δ features | Does trajectory matter? |
# | `full` | All combined | Upper bound |
#
# **Feature Types:** `cna`, `frag`, `motif`, `combined` (100 total combinations)

# %% [cell 34]
from sklearn.utils.class_weight import compute_sample_weight

df2=feat_df.dropna(subset=['CT']).copy()
df2['y']=(df2['CT']=='PR').astype(int)
y_s2=df2['y'].values; loocv=LeaveOneOut()
print(f'Stage 2: n={len(df2)} PR={y_s2.sum()} non-PR={(y_s2==0).sum()}')

# ── Engineer motif delta features (ใช้ใน change/full scenario) ───────────
for k in range(1, 11):
    feat_df[f'dmpc{k}_dt2']  = feat_df[f'mpc{k}_t2'] - feat_df[f'mpc{k}_t1']  # Δ T0→W3
    feat_df[f'dmpc{k}_dt3']  = feat_df[f'mpc{k}_t3'] - feat_df[f'mpc{k}_t1']  # Δ T0→W12
    feat_df[f'dmpc{k}_dt23'] = feat_df[f'mpc{k}_t3'] - feat_df[f'mpc{k}_t2']  # Δ W3→W12

# rebuild df2 หลัง engineer features ใหม่
df2 = feat_df.dropna(subset=['CT']).copy()
df2['y'] = (df2['CT']=='PR').astype(int)
y_s2 = df2['y'].values

def get_fc(df, scenario, ftype):
    """
    5 scenarios ใหม่ — แต่ละ scenario ใช้ข้อมูลเฉพาะช่วงเวลาของตัวเอง
    ไม่ overlap กัน (ยกเว้น full ที่รวมทั้งหมด)

    baseline : T0 เท่านั้น              → ทำนายได้ก่อนรักษาไหม?
    early    : W3-4 เท่านั้น            → สัญญาณช่วงแรกของการรักษา
    late     : W12 เท่านั้น             → สัญญาณตอนสิ้นสุดการรักษา
    change   : Δ(T0→W3), Δ(T0→W12),   → trajectory สำคัญแค่ไหน?
               Δ(W3→W12)
    full     : ทุกอย่างรวมกัน          → upper bound ของ model
    """
    c = []

    # ── CNA / ichorCNA features ──────────────────────────────────────────
    if ftype in ['cna', 'combined']:

        if scenario in ['baseline', 'full']:
            c += [f'{f}_t1' for f in ICHOR_FEATS]
            c += ['TF_pos_t1', 'cp_t1']

        if scenario in ['early', 'full']:
            c += [f'{f}_t2' for f in ICHOR_FEATS]
            c += ['TF_pos_t2', 'cp_t2']

        if scenario in ['late', 'full']:
            c += [f'{f}_t3' for f in ICHOR_FEATS]
            c += ['TF_pos_t3', 'cp_t3']

        if scenario in ['change', 'full']:
            c += [f'd_{f}_dt2'    for f in ICHOR_FEATS]   # Δ T0→W3
            c += [f'pct_{f}_dt2'  for f in ICHOR_FEATS]
            c += [f'd_{f}_dt3'    for f in ICHOR_FEATS]   # Δ T0→W12
            c += [f'pct_{f}_dt3'  for f in ICHOR_FEATS]
            c += [f'd_{f}_dt23'   for f in ICHOR_FEATS]   # Δ W3→W12
            c += [f'pct_{f}_dt23' for f in ICHOR_FEATS]
            c += ['TF_clear_t2', 'TF_clear_t3', 'TF_rebound_t3']
            c += ['dcp_dt2', 'dcp_dt3']

    # ── Fragment length features ──────────────────────────────────────────
    if ftype in ['frag', 'combined']:

        if scenario in ['baseline', 'full']:
            c += [f'{f}_t1' for f in FRAG_FEATS]

        if scenario in ['early', 'full']:
            c += [f'{f}_t2' for f in FRAG_FEATS]

        if scenario in ['late', 'full']:
            c += [f'{f}_t3' for f in FRAG_FEATS]

        if scenario in ['change', 'full']:
            c += [f'd_{f}_dt2'    for f in FRAG_FEATS]   # Δ T0→W3
            c += [f'pct_{f}_dt2'  for f in FRAG_FEATS]
            c += [f'd_{f}_dt3'    for f in FRAG_FEATS]   # Δ T0→W12
            c += [f'pct_{f}_dt3'  for f in FRAG_FEATS]
            c += [f'd_{f}_dt23'   for f in FRAG_FEATS]   # Δ W3→W12
            c += [f'pct_{f}_dt23' for f in FRAG_FEATS]

    # ── End-motif PCA features ────────────────────────────────────────────
    if ftype in ['motif', 'combined']:

        if scenario in ['baseline', 'full']:
            c += [f'mpc{k}_t1' for k in range(1, 11)]

        if scenario in ['early', 'full']:
            c += [f'mpc{k}_t2' for k in range(1, 11)]

        if scenario in ['late', 'full']:
            c += [f'mpc{k}_t3' for k in range(1, 11)]

        if scenario in ['change', 'full']:
            c += [f'dmpc{k}_dt2'  for k in range(1, 11)]  # Δ T0→W3
            c += [f'dmpc{k}_dt3'  for k in range(1, 11)]  # Δ T0→W12
            c += [f'dmpc{k}_dt23' for k in range(1, 11)]  # Δ W3→W12

    return [col for col in c if col in df.columns]

def make_s2():
    return {
        'LR': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                         ('clf',LogisticRegression(class_weight='balanced',max_iter=1000,random_state=SEED))]),
        'RF': Pipeline([('imp',SimpleImputer(strategy='median')),
                         ('clf',RandomForestClassifier(n_estimators=200,class_weight='balanced',random_state=42))]),
        'SVM':Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                         ('clf',SVC(probability=True,class_weight='balanced',random_state=SEED))]),
        'GB':   Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                          ('clf',GradientBoostingClassifier(n_estimators=50,random_state=42))]),
        'LGBM': Pipeline([('imp',SimpleImputer(strategy='median')),
                          ('clf',lgb.LGBMClassifier(
                              n_estimators=100,learning_rate=0.05,max_depth=4,
                              num_leaves=15,min_child_samples=3,
                              class_weight='balanced',random_state=42,
                              verbosity=-1,n_jobs=-1))]),
    }

# 5 scenarios × 4 feature sets × 5 models = 100 combinations
SCENARIOS    = ['baseline', 'early', 'late', 'change', 'full']
FEATURE_SETS = ['cna', 'frag', 'motif', 'combined']
records=[]; best_auc=0; best_cfg=None; sw_all=compute_sample_weight('balanced',y_s2)
set_global_seed(SEED)  # reset before scenario sweep

print(f"\n{'Mdl':4} {'Scenario':10} {'Feats':9} {'#F':>4} {'AUC':>6} {'Acc':>6} {'F1':>6}")
print('-'*58)
for sc in SCENARIOS:
    for ft in FEATURE_SETS:
        fc=get_fc(df2,sc,ft)
        if not fc: continue
        X=df2[fc].values
        for nm,mdl in make_s2().items():
            try:
                if nm in ('GB','LGBM'):
                    yp=cross_val_predict(mdl,X,y_s2,cv=loocv,method='predict_proba',
                                         params={'clf__sample_weight':sw_all})[:,1]
                else:
                    yp=cross_val_predict(mdl,X,y_s2,cv=loocv,method='predict_proba')[:,1]
                yd=(yp>0.5).astype(int)
                auc=roc_auc_score(y_s2,yp); acc=accuracy_score(y_s2,yd); f1v=f1_score(y_s2,yd,zero_division=0)
                rec={'model':nm,'scenario':sc,'feat_type':ft,'n_feat':len(fc),'AUC':auc,'Acc':acc,'F1':f1v,'yp':yp,'fc':fc}
                records.append(rec)
                if auc>best_auc: best_auc=auc; best_cfg=rec
                print(f"{nm:4} {sc:10} {ft:9} {len(fc):4d} {auc:6.3f} {acc:6.3f} {f1v:6.3f}")
            except Exception as e: print(f'  ERR {nm}/{sc}/{ft}: {e}')

s2_df=pd.DataFrame([{k:v for k,v in r.items() if k not in ['yp','fc']} for r in records])
y_true=y_s2
print(f"\nBest: {best_cfg['model']} | {best_cfg['scenario']} | {best_cfg['feat_type']} | AUC={best_auc:.3f}")

# %% [cell 35]
# ### Stage 2 · AUC Heatmap + Best Model Detail

# %% [cell 36]
# AUC heatmap — show all 5 models (including LGBM)
models_to_show = [m for m in ['LR', 'RF', 'SVM', 'GB', 'LGBM'] if m in s2_df['model'].unique()]
n_models = len(models_to_show)
fig, axes = plt.subplots(1, n_models, figsize=(5.5 * n_models, 5.5))
fig.patch.set_facecolor('white')
fig.suptitle('Stage 2 AUC: Scenario × Feature Set (LOOCV | PR vs non-PR)',
             fontsize=13, fontweight='bold', y=1.02)
if n_models == 1:
    axes = [axes]
for ax, nm in zip(axes, models_to_show):
    sub=s2_df[s2_df['model']==nm]
    pvt=sub.pivot_table(index='feat_type',columns='scenario',values='AUC')
    pvt=pvt.reindex(FEATURE_SETS).reindex(SCENARIOS,axis=1)
    sns.heatmap(pvt,annot=True,fmt='.3f',cmap='RdYlGn',ax=ax,
                vmin=0.4,vmax=1.0,linewidths=0.6,linecolor='white',annot_kws={'size':9,'weight':'bold'})
    ax.set_title(nm,fontweight='bold'); ax.set_xlabel('Scenario'); ax.tick_params(axis='x',rotation=30)
    ax.set_ylabel('Feature Set' if nm=='LR' else '')
add_panel_labels(axes)
plt.tight_layout(); plt.savefig(f'{OUT}/fig_s2_heatmap.png',dpi=150,bbox_inches='tight'); plt.show()

# Best model detail
fig,axes=plt.subplots(1,3,figsize=(18,5))
fig.suptitle(f"Best S2: {best_cfg['model']} | {best_cfg['scenario']} | {best_cfg['feat_type']} | AUC={best_auc:.3f}",
             fontsize=12,fontweight='bold')
fp2,tp2,_=roc_curve(y_true,best_cfg['yp'])
axes[0].plot(fp2,tp2,'#e74c3c',lw=2.5,label=f"AUC={best_auc:.3f}")
axes[0].plot([0,1],[0,1],'k--',lw=1); axes[0].set_xlabel('1-Specificity'); axes[0].set_ylabel('Sensitivity')
axes[0].set_title('ROC',fontweight='bold'); axes[0].legend(fontsize=10); axes[0].grid(alpha=0.3)
axes[0].spines[['top','right']].set_visible(False)

yd_b=(best_cfg['yp']>0.5).astype(int)
ConfusionMatrixDisplay(confusion_matrix(y_true,yd_b),display_labels=['SD+PD','PR']).plot(
    ax=axes[1],colorbar=False,cmap='Blues')
sen=recall_score(y_true,yd_b,zero_division=0); spe=recall_score(1-y_true,1-yd_b,zero_division=0)
axes[1].set_title(f"Confusion\nSens={sen:.2f} Spec={spe:.2f} F1={f1_score(y_true,yd_b,zero_division=0):.2f}",fontweight='bold')

bfc=best_cfg['fc']; Xfi=df2[bfc].fillna(df2[bfc].median())
rf_s2=RandomForestClassifier(n_estimators=300,class_weight='balanced',random_state=42).fit(Xfi,y_true)
imp_s2=pd.Series(rf_s2.feature_importances_,index=bfc).sort_values(ascending=False).head(15)
axes[2].barh(imp_s2.index[::-1],imp_s2.values[::-1],color='#e74c3c',alpha=0.85,edgecolor='white')
axes[2].set_title('Top-15 Feature Importance',fontweight='bold'); axes[2].set_xlabel('Importance')
axes[2].grid(axis='x',alpha=0.3); axes[2].spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig(f'{OUT}/fig_s2_detail.png',dpi=150,bbox_inches='tight'); plt.show()

# Calibration
fig,ax=plt.subplots(figsize=(5,5))
fp_,mp_=calibration_curve(y_true,best_cfg['yp'],n_bins=5)
ax.plot(mp_,fp_,'o-',lw=2,color='#e74c3c',ms=8,label='Model')
ax.plot([0,1],[0,1],'k--',lw=1.5,label='Perfect')
brier=brier_score_loss(y_true,best_cfg['yp'])
ax.set_title(f"Calibration (Brier={brier:.3f})",fontweight='bold')
ax.set_xlabel('Mean Predicted Prob'); ax.set_ylabel('Fraction Positives')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(f'{OUT}/fig_s2_calibration.png',dpi=150,bbox_inches='tight'); plt.show()
print(f'Brier={brier:.3f}  (0=perfect, 0.25=no-skill)')
feat_df.to_csv(f'{OUT}/feature_matrix.csv',index=False)
s2_df.to_csv(f'{OUT}/table_s2_results.csv',index=False)
print('Saved: feature_matrix.csv, table_s2_results.csv')

# %% [cell 37]
# ### Stage 2 · SHAP Feature Explanation (Best Model per Scenario)

# %% [cell 38]
# ── Stage 2 · SHAP Explanation (Improved / Cleaner Layout) ──────────────
import re
import shap
import lightgbm as lgb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC

# -----------------------------
# Inputs
# -----------------------------
best_fc2  = best_cfg['fc']
best_nm2  = best_cfg['model']
X2_raw    = df2[best_fc2].values
y2        = y_s2.copy()
n2        = len(X2_raw)

print(f'[SHAP S2] model={best_nm2}  features={len(best_fc2)}  n={n2}')

# -----------------------------

def normalize_binary_shap_output(sv):
    """
    Convert SHAP output to shape (n_samples, n_features)
    for binary classification.
    """
    if isinstance(sv, list):
        return sv[1]
    if hasattr(sv, 'ndim') and sv.ndim == 3:
        return sv[:, :, 1]
    return sv

def prettify_feature_name(name):
    """
    Shorten long feature names for cleaner plotting.
    """
    base_map = {
        'tumor_fraction'                        : 'tumor_fraction',
        'gc_mad'                                : 'gc_mad',
        'ploidy'                                : 'ploidy',
        'short_fraction_s150'                   : 'short_frac_s150',
        'short_to_long_100_150_over_151_220'    : 'S/L ratio',
        'frac_100_150'                          : 'frac_100_150',
        'frac_151_180'                          : 'frac_151_180',
        'frac_181_220'                          : 'frac_181_220',
        'frac_221_300'                          : 'frac_221_300',
        'mean'                                  : 'mean',
        'median'                                : 'median',
        'skewness'                              : 'skewness',
        'kurtosis_excess'                       : 'kurtosis',
        'TF_pos'                                : 'TF_pos',
        'TF_clear'                              : 'TF_clear',
        'cp'                                    : 'cp',
        'mpc1'                                  : 'mpc1',
        'mpc2'                                  : 'mpc2',
        'mpc3'                                  : 'mpc3',
        'mpc4'                                  : 'mpc4',
        'mpc5'                                  : 'mpc5',
        'mpc6'                                  : 'mpc6',
        'mpc7'                                  : 'mpc7',
        'mpc8'                                  : 'mpc8',
        'mpc9'                                  : 'mpc9',
    }

    # detect suffix เช่น _t1, _t2, _t3, _dt2, _dt3, _dt23
    m = re.match(r'^(.*?)(_(t1|t2|t3|dt2|dt3|dt23))$', name)
    if m:
        stem = m.group(1)
        suffix = m.group(2)
    else:
        stem = name
        suffix = ''

    pretty_stem = base_map.get(stem, stem)
    return pretty_stem + suffix

# -----------------------------
# Re-train best Stage-2 model on full data
# -----------------------------
imp2 = SimpleImputer(strategy='median')
sc2  = StandardScaler()

X2_imp = imp2.fit_transform(X2_raw)
sw2    = compute_sample_weight('balanced', y2)

if best_nm2 == 'LR':
    X2_sc = sc2.fit_transform(X2_imp)
    clf2  = LogisticRegression(class_weight='balanced', max_iter=1000)
    clf2.fit(X2_sc, y2)

    explainer2 = shap.LinearExplainer(
        clf2, X2_sc, feature_perturbation='interventional'
    )
    sv2     = explainer2.shap_values(X2_sc)
    shap2   = normalize_binary_shap_output(sv2)
    X2_plot = X2_sc

elif best_nm2 == 'RF':
    clf2 = RandomForestClassifier(
        n_estimators=200,
        class_weight='balanced',
        random_state=42
    )
    clf2.fit(X2_imp, y2)

    explainer2 = shap.TreeExplainer(clf2)
    sv2        = explainer2.shap_values(X2_imp)
    shap2      = normalize_binary_shap_output(sv2)
    X2_plot    = X2_imp

elif best_nm2 == 'GB':
    clf2 = GradientBoostingClassifier(n_estimators=50, random_state=42)
    clf2.fit(X2_imp, y2, sample_weight=sw2)

    explainer2 = shap.TreeExplainer(clf2)
    sv2        = explainer2.shap_values(X2_imp)
    shap2      = normalize_binary_shap_output(sv2)
    X2_plot    = X2_imp

elif best_nm2 == 'LGBM':
    clf2 = lgb.LGBMClassifier(
        n_estimators=100,
        learning_rate=0.05,
        max_depth=4,
        num_leaves=15,
        min_child_samples=3,
        class_weight='balanced',
        random_state=42,
        verbosity=-1
    )
    clf2.fit(X2_imp, y2, sample_weight=sw2)

    explainer2 = shap.TreeExplainer(clf2)
    sv2        = explainer2.shap_values(X2_imp)
    shap2      = normalize_binary_shap_output(sv2)
    X2_plot    = X2_imp

else:  # SVM — KernelExplainer
    X2_sc = sc2.fit_transform(X2_imp)
    clf2  = SVC(probability=True, class_weight='balanced', random_state=SEED)
    clf2.fit(X2_sc, y2)

    bg = shap.kmeans(X2_sc, min(20, n2))
    explainer2 = shap.KernelExplainer(clf2.predict_proba, bg)
    sv2        = explainer2.shap_values(X2_sc, nsamples=100)
    shap2      = normalize_binary_shap_output(sv2)
    X2_plot    = X2_sc

print(f'X2_plot shape: {X2_plot.shape}')
print(f'shap2 shape: {shap2.shape}')

assert shap2.shape[0] == X2_plot.shape[0], "sample dimension mismatch"
assert shap2.shape[1] == X2_plot.shape[1], "feature dimension mismatch"

# -----------------------------
# Prepare top features
# -----------------------------
mean_abs2 = np.abs(shap2).mean(axis=0)
top_n2    = min(15, len(best_fc2))
top_idx2  = np.argsort(mean_abs2)[-top_n2:][::-1]

top_feat2   = [best_fc2[i] for i in top_idx2]
top_label2  = [prettify_feature_name(f) for f in top_feat2]
top_shap2   = shap2[:, top_idx2]
top_X2      = X2_plot[:, top_idx2]
top_mean2   = mean_abs2[top_idx2]

shap_cmap = LinearSegmentedColormap.from_list(
    'shap_blue_pink',
    ['#1E88E5', '#FF0052']
)

# =============================
# Figure 1: Beeswarm + Bar Importance
# =============================
fig = plt.figure(figsize=(13.8, 7.2))
gs = fig.add_gridspec(
    1, 3,
    width_ratios=[2.15, 0.12, 1.70],
    wspace=0.65
)

ax_bee = fig.add_subplot(gs[0, 0])
cax    = fig.add_subplot(gs[0, 1])
ax_bar = fig.add_subplot(gs[0, 2])

fig.suptitle(
    f'Stage 2 SHAP — {best_nm2} | {best_cfg["scenario"]} | {best_cfg["feat_type"]}\n'
    f'(AUC={best_auc:.3f}  PR vs non-PR)',
    fontsize=15,
    fontweight='bold',
    y=1.00
)

# ---- Left: custom beeswarm
rng = np.random.default_rng(42)

for row in range(top_n2):
    sv = top_shap2[:, row]
    fv = top_X2[:, row]

    lo = np.nanpercentile(fv, 5)
    hi = np.nanpercentile(fv, 95)

    if hi > lo:
        fv_norm = np.clip((fv - lo) / (hi - lo), 0, 1)
    else:
        fv_norm = np.full_like(fv, 0.5, dtype=float)

    y = row + rng.uniform(-0.18, 0.18, size=len(sv))

    ax_bee.scatter(
        sv, y,
        c=fv_norm,
        cmap=shap_cmap,
        s=18,
        alpha=0.95,
        linewidths=0,
        zorder=3
    )

ax_bee.axvline(0, color='#888888', lw=1.2, ls='--', zorder=1)
ax_bee.set_yticks(range(top_n2))
ax_bee.set_yticklabels(top_label2, fontsize=10.5)
ax_bee.invert_yaxis()

ax_bee.set_title('Beeswarm (Top-15)', fontsize=13.5, fontweight='bold', pad=8)
ax_bee.set_xlabel('SHAP value (impact on model output)', fontsize=11.5)
ax_bee.grid(axis='y', alpha=0.18, linestyle=':', linewidth=0.8)
ax_bee.grid(axis='x', alpha=0.15, linestyle='--', linewidth=0.8)
ax_bee.spines[['top', 'right', 'left']].set_visible(False)
ax_bee.tick_params(axis='y', pad=6)
ax_bee.tick_params(axis='x', labelsize=10.5)

# ---- Colorbar
sm = plt.cm.ScalarMappable(cmap=shap_cmap)
sm.set_array([])

cbar = fig.colorbar(sm, cax=cax)
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['Low', 'High'])
cbar.set_label('')
cbar.ax.set_title('Feature value', fontsize=11, pad=8)
cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=10)

# ---- Right: mean |SHAP| importance
ypos = np.arange(top_n2)

ax_bar.barh(
    ypos,
    top_mean2,
    color='#1E88E5',
    height=0.68
)

ax_bar.set_yticks(ypos)
ax_bar.set_yticklabels(top_label2, fontsize=10.5)
ax_bar.invert_yaxis()

ax_bar.set_title('Mean |SHAP| Importance', fontsize=13.5, fontweight='bold', pad=8)
ax_bar.set_xlabel('mean |SHAP value|', fontsize=11.5)
ax_bar.grid(axis='x', alpha=0.18, linestyle='--', linewidth=0.8)
ax_bar.spines[['top', 'right', 'left']].set_visible(False)
ax_bar.tick_params(axis='y', pad=4)
ax_bar.tick_params(axis='x', labelsize=10.5)

fig.subplots_adjust(left=0.23, right=0.98, top=0.88, bottom=0.12)

add_panel_labels(axes)
plt.savefig(f'{OUT}fig_shap_stage2.png', dpi=150, bbox_inches='tight')
plt.show()
print('✓ fig_shap_stage2.png')

# -----------------------------
# Importance table
# -----------------------------
shap_imp_s2 = (
    pd.Series(np.abs(shap2).mean(axis=0), index=best_fc2)
    .sort_values(ascending=False)
)

print('\nTop-10 SHAP features (Stage 2):')
print(shap_imp_s2.head(10).to_string())

# =============================
# Figure 2: Directional Contribution
# =============================
pr_idx  = np.where(y2 == 1)[0]
npr_idx = np.where(y2 == 0)[0]

mean_pr2  = shap2[pr_idx].mean(axis=0)
mean_npr2 = shap2[npr_idx].mean(axis=0)

ord2 = np.argsort(np.abs(shap2).mean(axis=0))[::-1][:top_n2]
feat_ord2 = [prettify_feature_name(best_fc2[j]) for j in ord2]
y_pos2 = np.arange(len(ord2))

fig2, ax2 = plt.subplots(figsize=(10.2, max(5.2, len(ord2) * 0.55)))

ax2.barh(
    y_pos2 - 0.18,
    mean_pr2[ord2],
    height=0.34,
    color='#2F6FB0',
    alpha=0.88,
    label='PR',
    edgecolor='white'
)
ax2.barh(
    y_pos2 + 0.18,
    mean_npr2[ord2],
    height=0.34,
    color='#D97663',
    alpha=0.88,
    label='non-PR',
    edgecolor='white'
)

# symmetric x-axis for cleaner visual balance
xmax = np.max(np.abs(np.r_[mean_pr2[ord2], mean_npr2[ord2]]))
xmax = max(xmax * 1.15, 0.05)

ax2.axvline(0, color='black', lw=1.0)
ax2.set_xlim(-xmax, xmax)
ax2.set_yticks(y_pos2)
ax2.set_yticklabels(feat_ord2, fontsize=10.5)
ax2.invert_yaxis()

ax2.set_xlabel('Mean SHAP value  (+ → PR,  − → non-PR)', fontsize=11.5)
ax2.set_title('SHAP Direction: PR vs non-PR — Top 15', fontsize=14, fontweight='bold', pad=10)
ax2.legend(fontsize=10, loc='upper right', frameon=True)
ax2.spines[['top', 'right']].set_visible(False)
ax2.grid(axis='x', alpha=0.25, linestyle='--')

plt.tight_layout(rect=[0.08, 0.04, 0.98, 0.96], pad=1.5)
plt.savefig(f'{OUT}fig_shap_stage2_direction.png', dpi=150, bbox_inches='tight')
plt.show()
print('✓ fig_shap_stage2_direction.png')

# %% [cell 39]
# ---
# # Part 5 · Model Improvement Pipeline
#
# **Systematic optimization steps:**
#
# | Step | Method | Purpose |
# |------|--------|---------|
# | 5.1 | 3-Method Consensus (RF + MI + Mann-Whitney) | Rank features by discriminability |
# | 5.2 | Redundancy Removal (\|Spearman r\| > 0.90) | Reduce multicollinearity |
# | 5.3 | Top-K Sweep × 7 Models | Find optimal feature count |
# | 5.4 | Soft-Vote Ensemble (diverse top-3) | Aggregate predictions |
# | 5.5 | Summary Table | Compare all models |
# | 5.7 | Per-Patient Tracker | Individual prediction audit |
# | 5.8 | Kaplan-Meier PFS Analysis | Clinical validation |
# | 5.9 | SHAP Waterfall per Patient | Explainability |

# %% [cell 40]
df2=feat_df.dropna(subset=['CT']).copy()
df2['y']=(df2['CT']=='PR').astype(int); y=df2['y'].values
ALL_FC=[c for c in feat_df.columns if c not in ['pid','CT'] and '_z_' not in c and not c.startswith('rate_')]
Xall_imp=SimpleImputer(strategy='median').fit_transform(df2[ALL_FC].values)
Xall=pd.DataFrame(Xall_imp,columns=ALL_FC)
print(f'Samples={len(df2)} PR={y.sum()} non-PR={(y==0).sum()} Features={len(ALL_FC)}')

# %% [cell 41]
# ### Step 5.1 · 3-Method Consensus Feature Importance

# %% [cell 42]
print('[1/3] RF importance ...')
rf500=RandomForestClassifier(n_estimators=500,class_weight='balanced',random_state=42,n_jobs=-1)
rf500.fit(Xall,y); imp_rf=pd.Series(rf500.feature_importances_,index=ALL_FC)

print('[2/3] Mutual Information ...')
mi_vals=mutual_info_classif(Xall.values,y,random_state=42); imp_mi=pd.Series(mi_vals,index=ALL_FC)

print('[3/3] Mann-Whitney ...')
mw_rows=[]
for col in ALL_FC:
    g1=df2[df2['y']==1][col].dropna(); g2=df2[df2['y']==0][col].dropna()
    if len(g1)>1 and len(g2)>1:
        _,p=mannwhitneyu(g1,g2,alternative='two-sided')
        mw_rows.append({'feature':col,'p_value':p,'diff':g1.median()-g2.median()})
mw_df2=pd.DataFrame(mw_rows).sort_values('p_value')

rank_rf=imp_rf.rank(ascending=False); rank_mi=imp_mi.rank(ascending=False)
rank_mw=(mw_df2.set_index('feature')['p_value'].reindex(ALL_FC).rank(ascending=True).fillna(len(ALL_FC)))
avg_rank=(rank_rf+rank_mi+rank_mw)/3
top20_rank=avg_rank.sort_values().head(20); top20_cols=top20_rank.index.tolist()

lmap={'kurtosis_excess':'Kurtosis','skewness':'Skewness','short_to_long_100_150_over_151_220':'S/L Ratio',
      'tumor_fraction':'TF','ploidy':'Ploidy','gc_mad':'GC-MAD','mean':'FragMean','median':'FragMedian',
      'short_fraction_s150':'ShortFrac'}
def cl(c):
    tp_={'_t1':' @BL','_t2':' @W3','_t3':' @W12'}
    if c.startswith('dcp'): return 'ΔCancerProb'
    if c.startswith('cp_'):
        for k,v in tp_.items():
            if c.endswith(k): return f'CancerProb{v}'
    if c.startswith('mpc'):
        s=c
        for k,v in tp_.items(): s=s.replace(k,v)
        return s.replace('mpc','MotifPC')
    if 'TF_' in c:
        s=c.replace('TF_pos','TF+').replace('TF_clear','TFClr').replace('TF_rebound','TFReb')
        for k,v in tp_.items(): s=s.replace(k,v)
        return s
    pfx=''; c2=c
    if c.startswith('d_'):    pfx='D ';  c2=c[2:]
    elif c.startswith('pct_'): pfx='%D '; c2=c[4:]
    for k,v in tp_.items():
        if c2.endswith(k): b=c2[:-(len(k))]; return f'{pfx}{lmap.get(b,b.replace("_"," "))}{v}'
    return pfx+lmap.get(c2,c2.replace('_',' '))

print(f"\n{'#':>3}  {'Feature':45}  {'RF':>6}  {'MI':>6}  {'MW-p':>8}")
print('  '+'-'*70)
for i,(fn,_) in enumerate(top20_rank.items(),1):
    pv=mw_df2.set_index('feature')['p_value'].get(fn,1.0)
    print(f'  {i:3d}  {fn:45}  {imp_rf[fn]:6.4f}  {imp_mi[fn]:6.4f}  {pv:8.4f}')

# %% [cell 43]
# ### Step 5.2 · Remove Redundant Features |Spearman r| > 0.90

# %% [cell 44]
Xtop=df2[top20_cols].fillna(df2[top20_cols].median())
corr_mat=Xtop.corr(method='spearman').abs(); to_drop=set()
for i in range(len(top20_cols)):
    if top20_cols[i] in to_drop: continue
    for j in range(i+1,len(top20_cols)):
        if top20_cols[j] in to_drop: continue
        if corr_mat.iloc[i,j]>0.90:
            worse=top20_cols[i] if avg_rank[top20_cols[i]]>avg_rank[top20_cols[j]] else top20_cols[j]
            kept=top20_cols[i if worse==top20_cols[j] else j]
            to_drop.add(worse); print(f"  Drop '{worse}'  r={corr_mat.iloc[i,j]:.2f} with '{kept}'")
clean_cols=[c for c in top20_cols if c not in to_drop]
print(f'After dedup: {len(clean_cols)} features (removed {len(to_drop)})')
for i,c in enumerate(clean_cols,1): print(f'  {i:2d}. {cl(c):30}  [{c}]')

# %% [cell 45]
# ### Step 5.3 · Top-K Sweep × 6 Models + Ensemble

# %% [cell 46]
def make_sw():
    return {
        'LR C=0.1': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                               ('clf',LogisticRegression(C=0.1,class_weight='balanced',max_iter=2000,random_state=SEED))]),
        'LR C=0.5': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                               ('clf',LogisticRegression(C=0.5,class_weight='balanced',max_iter=2000,random_state=SEED))]),
        'RF d=4':   Pipeline([('imp',SimpleImputer(strategy='median')),
                               ('clf',RandomForestClassifier(n_estimators=300,max_depth=4,class_weight='balanced',random_state=42))]),
        'RF d=6':   Pipeline([('imp',SimpleImputer(strategy='median')),
                               ('clf',RandomForestClassifier(n_estimators=300,max_depth=6,class_weight='balanced',random_state=42))]),
        'GB d=2':   Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                               ('clf',GradientBoostingClassifier(n_estimators=60,max_depth=2,learning_rate=0.05,random_state=42))]),
        'SVM rbf':  Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                               ('clf',SVC(C=1.0,kernel='rbf',probability=True,class_weight='balanced',random_state=SEED))]),
        'LGBM':     Pipeline([('imp',SimpleImputer(strategy='median')),
                               ('clf',lgb.LGBMClassifier(
                                   n_estimators=150,learning_rate=0.05,max_depth=4,
                                   num_leaves=15,min_child_samples=3,
                                   class_weight='balanced',random_state=42,
                                   verbosity=-1,n_jobs=-1))]),
    }

def run_loocv_weighted(mdl, X, y, sw):
    loo = LeaveOneOut()
    yp  = np.zeros(len(y))
    for tr, te in loo.split(X):
        m = clone(mdl)
        m.fit(X[tr], y[tr], clf__sample_weight=sw[tr])
        yp[te[0]] = m.predict_proba(X[te])[:, 1]
    return yp

from sklearn.base import clone

k_list = [3, 5, 8, 10, 12, 15, len(clean_cols)]
sweep  = []
sw_all = compute_sample_weight('balanced', y)
set_global_seed(SEED)  # reset before sweep

print(f"  {'k':>3}  {'Model':12}  {'AUC':>6}  {'Acc':>6}  {'Sens':>6}  {'Spec':>6}  {'F1':>6}")
print('  ' + '-'*55)

for k in k_list:
    ck = clean_cols[:k]
    Xk = df2[ck].values
    for nm, mdl in make_sw().items():
        try:
            if nm in ('GB d=2', 'LGBM'):
                yp = run_loocv_weighted(mdl, Xk, y, sw_all)
            else:
                yp = cross_val_predict(mdl, Xk, y, cv=LeaveOneOut(),
                                       method='predict_proba')[:, 1]
            yd  = (yp > 0.5).astype(int)
            auc = roc_auc_score(y, yp)
            acc = accuracy_score(y, yd)
            f1v = f1_score(y, yd, zero_division=0)
            sen = recall_score(y, yd, zero_division=0)
            spe = recall_score(1-y, 1-yd, zero_division=0)
            sweep.append({'k':k,'model':nm,'AUC':auc,'Acc':acc,'F1':f1v,
                          'Sens':sen,'Spec':spe,'yp':yp,'cols':ck})
            print(f'  {k:3d}  {nm:12}  {auc:6.3f}  {acc:6.3f}  {sen:6.3f}  {spe:6.3f}  {f1v:6.3f}')
        except Exception as e:
            print(f'  ERR k={k} {nm}: {e}')

best_single = max(sweep, key=lambda r: r['AUC'])
print(f"\nBest single: k={best_single['k']}  {best_single['model']}  AUC={best_single['AUC']:.3f}")

# AUC heatmap — improved styling
swdf  = pd.DataFrame([{k2:v for k2,v in r.items() if k2 not in ['yp','cols']} for r in sweep])
pivot = swdf.pivot_table(index='model', columns='k', values='AUC')

fig, ax = plt.subplots(figsize=(14, 5.5))
fig.patch.set_facecolor('white')
sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn', ax=ax, vmin=0.5, vmax=1.0,
            linewidths=0.8, linecolor='white',
            annot_kws={'size': 10, 'weight': 'bold'},
            cbar_kws={'label': 'AUC (LOOCV)', 'shrink': 0.85})
ax.set_title('AUC — Model × Top-K Features  (LOOCV | PR vs Non-PR)',
             fontsize=13, fontweight='bold', pad=12)
ax.set_xlabel('K features (top by consensus rank)', fontsize=10)
ax.set_ylabel('Model', fontsize=10)
ax.tick_params(axis='x', rotation=0)
if best_single['model'] in pivot.index and best_single['k'] in pivot.columns:
    ri = list(pivot.index).index(best_single['model'])
    ci = list(pivot.columns).index(best_single['k'])
    ax.add_patch(plt.Rectangle((ci,ri),1,1,fill=False,edgecolor='black',lw=3,zorder=5))
    ax.text(ci+0.5, ri+0.5, '★', ha='center', va='center',
            fontsize=16, fontweight='bold', zorder=6)
add_panel_labels(ax)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/figB_sweep_heatmap.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [cell 47]
# ### Step 5.4 · Diverse Soft-Vote Ensemble + Final Visualisation

# %% [cell 48]
best_k=best_single['k']; best_cols=clean_cols[:best_k]
def pdiv3(res,k):
    at_k=[r for r in res if r['k']==k]; fams={'LR':[],'RF':[],'GB':[],'SVM':[],'LGBM':[]}
    for r in at_k:
        for f in fams:
            if r['model'].startswith(f) or r['model']==f: fams[f].append(r); break
    best_f=[max(v,key=lambda r:r['AUC']) for v in fams.values() if v]
    return sorted(best_f,key=lambda r:-r['AUC'])[:3]

top3=pdiv3(sweep,best_k); top3_names=[r['model'] for r in top3]
print(f'Ensemble (k={best_k}):')
for r in top3: print(f"  {r['model']:12}  AUC={r['AUC']:.3f}")

sw_e=compute_sample_weight('balanced',y); ens_p=np.zeros(len(y)); Xb=df2[best_cols].values
set_global_seed(SEED)  # reset before ensemble LOOCV
for tr,te in LeaveOneOut().split(Xb):
    Xtr_,Xte_=Xb[tr],Xb[te]; ytr_=y[tr]; sw_tr=sw_e[tr]; probs=[]
    for nm in top3_names:
        m=make_sw()[nm]
        if nm in ('GB d=2','LGBM'): m.fit(Xtr_,ytr_,**{'clf__sample_weight':sw_tr})
        else:                       m.fit(Xtr_,ytr_)
        probs.append(m.predict_proba(Xte_)[0,1])
    ens_p[te[0]]=np.mean(probs)

ens_auc=roc_auc_score(y,ens_p); ens_yd=(ens_p>0.5).astype(int)
ens_acc=accuracy_score(y,ens_yd); ens_f1=f1_score(y,ens_yd,zero_division=0)
ens_sens=recall_score(y,ens_yd,zero_division=0); ens_spec=recall_score(1-y,1-ens_yd,zero_division=0)
print(f'Ensemble AUC={ens_auc:.3f} Acc={ens_acc:.3f} Sens={ens_sens:.3f} Spec={ens_spec:.3f} F1={ens_f1:.3f}')

orig_cols=[c for c in ['tumor_fraction_t1','ploidy_t1','gc_mad_t1','TF_pos_t1',
                        'd_tumor_fraction_dt2','d_gc_mad_dt2','TF_pos_t2','TF_clear_t2'] if c in df2.columns]
p_orig=Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                  ('clf',RandomForestClassifier(n_estimators=300,class_weight='balanced',random_state=42))])
yp_orig=cross_val_predict(p_orig,df2[orig_cols].values,y,cv=LeaveOneOut(),method='predict_proba')[:,1]
auc_orig=roc_auc_score(y,yp_orig)

best_yp=ens_p if ens_auc>=best_single['AUC'] else best_single['yp']
best_auc_f=max(ens_auc,best_single['AUC'])
best_lbl=f'Ensemble(k={best_k})' if ens_auc>=best_single['AUC'] else f"{best_single['model']}(k={best_k})"
print(f'Original RF: AUC={auc_orig:.3f}  Best: AUC={best_auc_f:.3f}  ΔAUC=+{best_auc_f-auc_orig:.3f}')

# ROC comparison
pal7=['#e74c3c','#3498db','#2ecc71','#f39c12','#9b59b6','#1abc9c','#7f8c8d']
fig,axes=plt.subplots(1,2,figsize=(13,5))
fig.suptitle('ROC Curves: Baseline Model versus Optimised Model Comparison', fontsize=13, fontweight='bold')
ax=axes[0]
for rec,col in zip(sorted([r for r in sweep if r['k']==best_k],key=lambda r:-r['AUC']),pal7):
    fp_,tp_,_=roc_curve(y,rec['yp']); ax.plot(fp_,tp_,color=col,lw=1.8,label=f"{rec['model']} {rec['AUC']:.3f}")
ax.plot([0,1],[0,1],'k--',lw=1); ax.set_xlabel('1 − Specificity  (False Positive Rate)'); ax.set_ylabel('Sensitivity  (True Positive Rate)')
ax.set_title(f'All Models  (k = {best_k} features)',fontweight='bold'); ax.legend(loc='lower right',fontsize=8); ax.grid(alpha=0.3)
ax.spines[['top','right']].set_visible(False)
ax=axes[1]
fpo,tpo,_=roc_curve(y,yp_orig); fps,tps,_=roc_curve(y,best_single['yp']); fpe,tpe,_=roc_curve(y,ens_p)
ax.plot(fpo,tpo,'#95a5a6',lw=2,ls='--',label=f'Original {auc_orig:.3f}')
ax.plot(fps,tps,'#3498db',lw=2,label=f"{best_single['model']} {best_single['AUC']:.3f}")
ax.plot(fpe,tpe,'#e74c3c',lw=2.5,label=f'Ensemble {ens_auc:.3f}')
ax.fill_between(fpe,tpe,alpha=0.07,color='#e74c3c')
ax.plot([0,1],[0,1],'k--',lw=1); ax.set_xlabel('1 − Specificity  (False Positive Rate)'); ax.set_ylabel('Sensitivity  (True Positive Rate)')
ax.set_title('Baseline Model versus Optimised Model',fontweight='bold'); ax.legend(loc='lower right',fontsize=9); ax.grid(alpha=0.3)
ax.annotate(f'+{best_auc_f-auc_orig:.3f} AUC',xy=(0.35,0.65),fontsize=13,color='#e74c3c',fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3',facecolor='#ffeaa7',alpha=0.85))
ax.spines[['top','right']].set_visible(False)
add_panel_labels(axes)
plt.tight_layout(); plt.savefig(f'{OUT}/figC_ROC.png',dpi=150,bbox_inches='tight'); plt.show()

# Per-patient scatter
rng2=np.random.default_rng(0); yd_b=(best_yp>0.5).astype(int); df2['correct']=(yd_b==y)
fig,axes=plt.subplots(1,2,figsize=(13,5))
fig.suptitle(f'Predicted Probability of Partial Response — {best_lbl}  (AUC = {best_auc_f:.3f})', fontsize=13, fontweight='bold')
ax=axes[0]
for i,ct in enumerate(['PR','SD','PD']):
    mask=df2['CT'].values==ct; prb=best_yp[mask]; jit=rng2.uniform(-0.14,0.14,len(prb))
    ax.scatter([i+1+j for j in jit],prb,color=C[ct],s=90,alpha=0.75,edgecolors='white',lw=0.8,zorder=3)
    ax.plot([i+0.78,i+1.22],[np.median(prb),np.median(prb)],color=C[ct],lw=3.5,solid_capstyle='round',zorder=4)
ax.axhline(0.5,color='black',ls='--',lw=1.5); ax.set_xticks([1,2,3])
n1=(df2['CT']=='PR').sum(); n2=(df2['CT']=='SD').sum(); n3=(df2['CT']=='PD').sum()
ax.set_xticklabels([f'PR(n={n1})',f'SD(n={n2})',f'PD(n={n3})']); ax.set_ylabel('Predicted P(PR)')
ax.set_ylim(-0.05,1.05); ax.set_title('Predicted P(PR) per Patient  (line = median)',fontweight='bold')
ax.grid(axis='y',alpha=0.3); ax.spines[['top','right']].set_visible(False)
ax=axes[1]
for i,ct in enumerate(['PR','SD','PD']):
    mask=df2['CT'].values==ct
    for prb_,cor_ in zip(best_yp[mask],df2['correct'].values[mask]):
        ax.scatter(i+1,prb_,marker='o' if cor_ else 'X',color=C[ct] if cor_ else '#e74c3c',
                   s=120,alpha=0.85,edgecolors='white',lw=0.8,zorder=3)
ax.axhline(0.5,color='black',ls='--',lw=1.5); ax.set_xticks([1,2,3]); ax.set_xticklabels(['PR','SD','PD'])
ax.set_ylabel('Predicted P(PR)'); ax.set_ylim(-0.05,1.05); ax.set_title('Classification Outcome: Correct (○) versus Misclassified (×)',fontweight='bold')
nc_=df2['correct'].sum(); ax.text(0.97,0.04,f'Acc={nc_}/{len(df2)}({nc_/len(df2)*100:.0f}%)',
    transform=ax.transAxes,ha='right',fontsize=11,fontweight='bold',bbox=dict(boxstyle='round',facecolor='#dfe6e9',alpha=0.9))
ax.grid(axis='y',alpha=0.3); ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig(f'{OUT}/figD_scatter.png',dpi=150,bbox_inches='tight'); plt.show()

# %% [cell 49]
# ### Step 5.5 · Summary Table & Final Results

# %% [cell 50]
rows_=[]
yd_o=(yp_orig>0.5).astype(int)
rows_.append({'Model':'Original RF (baseline)','k':len(orig_cols),'AUC':round(auc_orig,3),
              'Acc':round(accuracy_score(y,yd_o),3),'Sens':round(recall_score(y,yd_o,zero_division=0),3),
              'Spec':round(recall_score(1-y,1-yd_o,zero_division=0),3),'F1':round(f1_score(y,yd_o,zero_division=0),3),'Note':'Baseline'})
for r in sorted(sweep,key=lambda r:-r['AUC'])[:6]:
    yd_=(r['yp']>0.5).astype(int)
    rows_.append({'Model':r['model'],'k':r['k'],'AUC':round(r['AUC'],3),
                  'Acc':round(accuracy_score(y,yd_),3),'Sens':round(recall_score(y,yd_,zero_division=0),3),
                  'Spec':round(recall_score(1-y,1-yd_,zero_division=0),3),'F1':round(f1_score(y,yd_,zero_division=0),3),'Note':'Single'})
rows_.append({'Model':f'Ensemble top-3 k={best_k}','k':best_k,'AUC':round(ens_auc,3),
              'Acc':round(ens_acc,3),'Sens':round(ens_sens,3),'Spec':round(ens_spec,3),'F1':round(ens_f1,3),'Note':'Ensemble'})

sum_df=pd.DataFrame(rows_).sort_values('AUC',ascending=False)
sum_df.to_csv(f'{OUT}/table_final_results.csv',index=False)
display(sum_df.style
        .background_gradient(subset=['AUC','Acc','F1'],cmap='RdYlGn',vmin=0.4,vmax=1.0)
        .format({'AUC':'{:.3f}','Acc':'{:.3f}','Sens':'{:.3f}','Spec':'{:.3f}','F1':'{:.3f}'})
        .set_caption('Model Comparison — LOOCV | PR vs Non-PR'))

print('='*55)
print('PIPELINE COMPLETE')
print('='*55)
print(f'  Stage 1 Best : {best_s1}  AUC={s1_res[best_s1]["AUC"]:.3f}')
print(f'  Original RF  : AUC={auc_orig:.3f}')
print(f'  Best improved: AUC={best_auc_f:.3f}  [{best_lbl}]')
print(f'  ΔAUC         : +{best_auc_f-auc_orig:.3f}')
print(f'  Final features: {clean_cols}')

# %% [cell 51]
# ### Step 5.6 · Learning Curve Analysis
#
# **Objective:** Diagnose bias-variance trade-off — is the model limited by data size or model complexity?
#
# - 🔵 Training AUC should approach 1.0 for complex models
# - 🔴 Validation AUC shows true generalization
# - Large gap → overfitting → consider regularization
# - Converging curves → underfitting or good fit

# %% [cell 52]
# ── Step 5.6 · Learning Curve Analysis (Best Model) ────────────────────────

from sklearn.model_selection import learning_curve

best_pipe_lc = make_sw().get(best_single['model'],
                              make_sw()['LR C=0.5'])  # fallback
Xb_lc = df2[best_cols].values

train_sizes_abs, train_scores, val_scores = learning_curve(
    best_pipe_lc, Xb_lc, y,
    cv=5, scoring='roc_auc',
    train_sizes=np.linspace(0.3, 1.0, 8),
    n_jobs=-1
)

train_mean = train_scores.mean(axis=1)
train_std  = train_scores.std(axis=1)
val_mean   = val_scores.mean(axis=1)
val_std    = val_scores.std(axis=1)

fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_facecolor('white')

ax.plot(train_sizes_abs, train_mean, 'o-', color='#3498db',
        lw=2, label='Training AUC', zorder=3)
ax.fill_between(train_sizes_abs,
                train_mean - train_std,
                train_mean + train_std,
                alpha=0.15, color='#3498db')

ax.plot(train_sizes_abs, val_mean, 'o-', color='#e74c3c',
        lw=2, label='Validation AUC (CV)', zorder=3)
ax.fill_between(train_sizes_abs,
                val_mean - val_std,
                val_mean + val_std,
                alpha=0.15, color='#e74c3c')

ax.axhline(0.5, color='#888888', ls=':', lw=1, label='Random baseline')
ax.set_xlabel('Training set size (samples)', fontsize=11)
ax.set_ylabel('AUC (ROC)', fontsize=11)
ax.set_title(f'Step 5.6 · Learning Curve — {best_single["model"]} (k={best_k} features)\n'
             f'Shaded area = ±1 SD across CV folds',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
ax.set_ylim(0.3, 1.05)
ax.grid(alpha=0.3)
ax.spines[['top','right']].set_visible(False)

gap = train_mean[-1] - val_mean[-1]
ax.text(0.98, 0.08,
        f'Final gap: {gap:+.3f}\n({"Overfitting" if gap > 0.15 else "Good fit" if gap < 0.05 else "Slight overfit"})',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=9, fontweight='bold',
        color='#e74c3c' if gap > 0.15 else '#27ae60',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc', alpha=0.9))

add_panel_labels(ax)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/figLC_learning_curve.png', dpi=150, bbox_inches='tight')
plt.show()
print(f'Learning curve saved. Train-Val gap at full data: {gap:+.3f}')

# %% [cell 53]
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.ticker as mticker

# ── รวบรวม models ที่จะแสดง ────────────────────────────────────────────────
# 1. Original RF baseline
yd_o   = (yp_orig > 0.5).astype(int)
# 2. Best single model
yd_bs  = (best_single['yp'] > 0.5).astype(int)
# 3. Ensemble
yd_ens = (ens_p > 0.5).astype(int)

cm_entries = [
    ('Original RF\n(Baseline)',   yp_orig,          yd_o,   '#95a5a6'),
    (f"{best_single['model']}\n(Best Single k={best_single['k']})",
                                   best_single['yp'], yd_bs,  '#3498db'),
    (f'Ensemble top-3\n(k={best_k})',
                                   ens_p,             yd_ens, '#e74c3c'),
]

labels = ['non-PR', 'PR']

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(
    'Confusion Matrices: Final Model Results  (LOOCV | PR vs Non-PR)',
    fontsize=13, fontweight='bold', y=1.02
)

for ax, (title, yp_, yd_, col) in zip(axes, cm_entries):
    cm = confusion_matrix(y, yd_)
    disp = ConfusionMatrixDisplay(cm, display_labels=labels)
    disp.plot(ax=ax, colorbar=False, cmap='Blues')

    # ── metric annotations ────────────────────────────────────────
    tn, fp, fn, tp = cm.ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv  = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv  = tn / (tn + fn) if (tn + fn) > 0 else 0
    acc  = (tp + tn) / len(y)
    auc_ = roc_auc_score(y, yp_)
    f1_  = f1_score(y, yd_, zero_division=0)

    metric_txt = (
        f'AUC  = {auc_:.3f}\n'
        f'Acc  = {acc:.3f}\n'
        f'Sens = {sens:.3f}   PPV = {ppv:.3f}\n'
        f'Spec = {spec:.3f}   NPV = {npv:.3f}\n'
        f'F1   = {f1_:.3f}'
    )
    ax.text(
        0.5, -0.28, metric_txt,
        transform=ax.transAxes, ha='center', va='top',
        fontsize=9, family='monospace',
        bbox=dict(boxstyle='round,pad=0.45', facecolor='#f8f9fa',
                  edgecolor=col, linewidth=1.8)
    )

    # ── title with color accent ────────────────────────────────────
    ax.set_title(title, fontweight='bold', fontsize=10, color=col, pad=10)
    ax.set_xlabel('Predicted label', fontsize=9)
    ax.set_ylabel('True label', fontsize=9)
    ax.spines[['top','right']].set_visible(False)

add_panel_labels(axes)
plt.tight_layout(rect=[0, 0.12, 1, 1])
plt.savefig(f'{OUT}/figE_confusion_final.png', dpi=150, bbox_inches='tight')
plt.show()
print('\n── Confusion Matrix Summary ──')
print(f'  {"Model":<30}  {"TP":>4}  {"TN":>4}  {"FP":>4}  {"FN":>4}  {"Sens":>6}  {"Spec":>6}  {"PPV":>6}  {"NPV":>6}')
print('  ' + '-'*75)
for title, yp_, yd_, _ in cm_entries:
    cm  = confusion_matrix(y, yd_)
    tn, fp, fn, tp = cm.ravel()
    sens = tp/(tp+fn) if (tp+fn)>0 else 0
    spec = tn/(tn+fp) if (tn+fp)>0 else 0
    ppv  = tp/(tp+fp) if (tp+fp)>0 else 0
    npv  = tn/(tn+fn) if (tn+fn)>0 else 0
    lbl  = title.replace('\n',' ')
    print(f'  {lbl:<30}  {tp:>4}  {tn:>4}  {fp:>4}  {fn:>4}  {sens:>6.3f}  {spec:>6.3f}  {ppv:>6.3f}  {npv:>6.3f}')

# %% [cell 54]
# ### Step 5.7 · Per-Patient Prediction Tracker (LOOCV | PR vs SD+PD)
#
# Audit individual predictions: shows P(PR) per patient with correct/misclassified status.

# %% [cell 55]
# ── Step 5.7 · Per-Patient Prediction Tracker (Best Model Only) ────────────

yd_best = (best_yp > 0.5).astype(int)
label_map = {1: 'PR', 0: 'non-PR'}

tracker = df2[['pid', 'CT']].copy().reset_index(drop=True)
tracker['prob']   = best_yp.round(3)
tracker['pred']   = [label_map[v] for v in yd_best]
tracker['actual'] = ['PR' if ct == 'PR' else 'non-PR' for ct in tracker['CT']]
tracker['ok']     = (tracker['pred'] == tracker['actual'])
tracker['result'] = tracker['ok'].map({True: '✓ Correct', False: '✗ Wrong'})

ct_order = {'PR': 0, 'SD': 1, 'PD': 2}
tracker['_o'] = tracker['CT'].map(ct_order)
tracker = tracker.sort_values(['_o', 'pid']).drop(columns='_o').reset_index(drop=True)

# ── Print ──────────────────────────────────────────────────────────────────
print(f"Best model: {best_lbl}  (AUC={best_auc_f:.3f})")
print(f"{'─'*62}")
print(f"  {'PID':<12}  {'CT':>5}  {'P(PR)':>7}  {'Pred':>7}  {'Result'}")
print(f"{'─'*62}")
prev_ct = None
for _, row in tracker.iterrows():
    if row['CT'] != prev_ct:
        print(f"  ── {row['CT']} ──")
        prev_ct = row['CT']
    icon = '✓' if row['ok'] else '✗'
    print(f"  {str(row['pid']):<12}  {row['CT']:>5}  {row['prob']:>7.3f}  {row['pred']:>7}  {icon}")
print(f"{'─'*62}")
n_ok = tracker['ok'].sum()
print(f"  Accuracy: {n_ok}/{len(tracker)} ({n_ok/len(tracker)*100:.1f}%)")
n_mis = (~tracker['ok']).sum()
if n_mis > 0:
    print(f"  Misclassified ({n_mis} patients):")
    for _, r in tracker[~tracker['ok']].iterrows():
        print(f"    PID={r['pid']}  CT={r['CT']}  P(PR)={r['prob']:.3f}  → predicted {r['pred']}")

# ── Plot ───────────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

C_CT = {'PR': '#4C78A8', 'SD': '#F58518', 'PD': '#E45756'}
cts  = tracker['CT'].values
pids = tracker['pid'].astype(str).values
probs = tracker['prob'].values
ok    = tracker['ok'].values
n_pat = len(tracker)
xs    = np.arange(n_pat)

fig, ax = plt.subplots(figsize=(14, 5.5))
fig.patch.set_facecolor('#f8f9fa')

# bars — green=correct, red=wrong
bar_cols = ['#27ae60' if o else '#e74c3c' for o in ok]
ax.bar(xs, probs, color=bar_cols, alpha=0.82, width=0.65,
       edgecolor='white', linewidth=0.5, zorder=3)

# threshold line
ax.axhline(0.5, color='#1a1a2e', ls='--', lw=1.4, zorder=4, label='threshold = 0.5')

# CT color band strip at top
for xi, ct in enumerate(cts):
    ax.axvspan(xi-0.5, xi+0.5, ymin=0.93, ymax=1.0,
               color=C_CT[ct], alpha=0.9, zorder=5)

# CT group dividers + labels
prev = None
group_starts = []
for xi, ct in enumerate(cts):
    if ct != prev:
        if xi > 0:
            ax.axvline(xi-0.5, color='#444', lw=1.2, zorder=6)
        group_starts.append((xi, ct))
        prev = ct
for i, (xi_s, ct) in enumerate(group_starts):
    xi_e = group_starts[i+1][0]-1 if i+1 < len(group_starts) else n_pat-1
    ax.text((xi_s+xi_e)/2, 1.045, ct, ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=C_CT[ct],
            transform=ax.get_xaxis_transform())

# prob labels on bars
for xi, (prob, o) in enumerate(zip(probs, ok)):
    ax.text(xi, prob + 0.02, f'{prob:.2f}', ha='center', va='bottom',
            fontsize=7, color='#333', rotation=0)

# X ticks
ax.set_xticks(xs)
ax.set_xticklabels(pids, rotation=55, ha='right', fontsize=8)
ax.set_xlim(-0.7, n_pat - 0.3)
ax.set_ylim(0, 1.12)
ax.set_ylabel('P(PR)', fontsize=11)
ax.set_xlabel('Patient ID', fontsize=11)
ax.set_title(f'Per-Patient Prediction — {best_lbl}  (AUC={best_auc_f:.3f})',
             fontsize=12, fontweight='bold', pad=14)
ax.grid(axis='y', alpha=0.25, zorder=0)
ax.spines[['top', 'right']].set_visible(False)

# legends
p1 = mpatches.Patch(color='#27ae60', alpha=0.85, label=f'Correct ({n_ok})')
p2 = mpatches.Patch(color='#e74c3c', alpha=0.85, label=f'Misclassified ({n_mis})')
ct_patches = [mpatches.Patch(color=v, label=k) for k, v in C_CT.items()]
ax.legend(handles=[p1, p2] + ct_patches, fontsize=8.5,
          loc='upper right', framealpha=0.9)

add_panel_labels(ax)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}/figF_patient_tracker.png', dpi=160, bbox_inches='tight',
            facecolor='#f8f9fa')
plt.show()

# ── Styled DataFrame ───────────────────────────────────────────────────────
display(tracker[['pid','CT','prob','pred','actual','result']].style
    .applymap(lambda v: 'color:#27ae60;font-weight:bold' if 'Correct' in str(v)
              else ('color:#e74c3c;font-weight:bold' if 'Wrong' in str(v) else ''),
              subset=['result'])
    .background_gradient(subset=['prob'], cmap='RdYlGn', vmin=0, vmax=1)
    .set_caption(f'Per-Patient Prediction vs CT Outcome — {best_lbl}'))

# %% [cell 56]
# ── Step 5.9 · Per-Patient Feature Contribution (True SHAP) ───────────────────

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import TwoSlopeNorm
import matplotlib.gridspec as gridspec
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut

# ── 1. เตรียม data & ระบุ best model ──────────────────────────────────────
USE_COLS = best_cols          # features ที่ best model ใช้ (จาก cell 46)
X_raw    = df2[USE_COLS].values
pids     = df2['pid'].values
cts      = df2['CT'].values
n_pat    = len(df2)

# ── Feature display name mapping ──────────────────────────────────────────
def short_name(col):
    m = {
        'tumor_fraction': 'TF', 'gc_mad': 'GC-MAD', 'ploidy': 'Ploidy',
        'short_to_long_100_150_over_151_220': 'S/L ratio',
        'short_fraction_s150': 'Short frac',
        'frac_100_150': 'Frac 100-150', 'frac_151_180': 'Frac 151-180',
        'frac_181_220': 'Frac 181-220', 'frac_221_300': 'Frac 221-300',
        'skewness': 'Skewness', 'kurtosis_excess': 'Kurtosis',
        'mean': 'Mean len', 'median': 'Median len',
        'TF_pos': 'TF-pos', 'TF_clear': 'TF-clear',
    }
    # ดึง base name + timepoint suffix
    for base, short in m.items():
        if col.startswith(base):
            suf = col[len(base):]
            suf = suf.replace('_t1','@T1').replace('_t2','@T2').replace('_t3','@T3')
            suf = suf.replace('_dt2','Δ@T2').replace('_dt3','Δ@T3')
            suf = suf.replace('_','')
            return f'{short}{suf}'
    # fallback: truncate
    return col[:22]

feat_labels = [short_name(c) for c in USE_COLS]
n_feat      = len(USE_COLS)

# ── 2. ระบุ best model type ────────────────────────────────────────────────
is_ensemble = ens_auc >= best_single['AUC']
model_type  = 'ensemble' if is_ensemble else best_single['model']
print(f"Best model: {best_lbl}  (AUC={best_auc_f:.3f})")
print(f"Features used: {n_feat}  → {USE_COLS}")

# ── 3. คำนวณ LOOCV contribution per patient ────────────────────────────────

def build_model(name):
    if 'LR' in name:
        C_ = 0.1 if '0.1' in name else 0.5
        return Pipeline([('imp', SimpleImputer(strategy='median')),
                         ('sc',  StandardScaler()),
                         ('clf', LogisticRegression(C=C_, class_weight='balanced',
                                                    max_iter=2000))])
    if 'RF' in name:
        d_ = 4 if 'd=4' in name else 6
        return Pipeline([('imp', SimpleImputer(strategy='median')),
                         ('clf', RandomForestClassifier(n_estimators=300, max_depth=d_,
                                  class_weight='balanced', random_state=42))])
    if 'GB' in name:
        return Pipeline([('imp', SimpleImputer(strategy='median')),
                         ('sc',  StandardScaler()),
                         ('clf', GradientBoostingClassifier(n_estimators=60, max_depth=2,
                                  learning_rate=0.05, random_state=42))])
    if 'SVM' in name:
        return Pipeline([('imp', SimpleImputer(strategy='median')),
                         ('sc',  StandardScaler()),
                         ('clf', SVC(C=1.0, kernel='rbf', probability=True,
                                     class_weight='balanced'))])
    if 'LGBM' in name:
        return Pipeline([('imp', SimpleImputer(strategy='median')),
                         ('clf', lgb.LGBMClassifier(
                              n_estimators=150, learning_rate=0.05, max_depth=4,
                              num_leaves=15, min_child_samples=3,
                              class_weight='balanced', random_state=42,
                              verbosity=-1, n_jobs=-1))])
    return None

def get_contribution_linear(pipe, X_te):
    """LR: coef × x_scaled  (signed contribution in log-odds space)"""
    imp   = pipe.named_steps['imp']
    sc    = pipe.named_steps['sc']
    clf   = pipe.named_steps['clf']
    X_imp = imp.transform(X_te)
    X_sc  = sc.transform(X_imp)
    coef  = clf.coef_[0]          # shape (n_feat,)
    # contribution ใน log-odds; scale ~พอดีกับ prob
    contribs = coef * X_sc[0]     # per-feature
    return contribs

def get_contribution_shap(pipe, X_tr, X_te, model_nm):
    """True SHAP contribution per patient (LOOCV-safe).
    X_tr, X_te are raw (pre-imputation) arrays from the pipeline split."""
    imp_s = pipe.named_steps['imp']
    X_tr_imp = imp_s.transform(X_tr)
    X_te_imp = imp_s.transform(X_te)

    clf = pipe.named_steps['clf']

    if 'LR' in model_nm:
        sc_s = pipe.named_steps['sc']
        X_bg = sc_s.transform(X_tr_imp)
        X_te_sc = sc_s.transform(X_te_imp)
        exp = shap.LinearExplainer(clf, X_bg, feature_perturbation='interventional')
        sv  = exp.shap_values(X_te_sc)
        return sv[0] if sv.ndim == 1 else sv[0]

    elif any(k in model_nm for k in ('RF','GB','LGBM')):
        exp = shap.TreeExplainer(clf, data=X_tr_imp,
                                feature_perturbation='interventional')
        sv  = exp.shap_values(X_te_imp)
        if isinstance(sv, list):  sv = sv[1]
        elif sv.ndim == 3:        sv = sv[:,:,1]
        return sv[0]

    else:  # SVM — fallback to permutation
        return get_contribution_permutation_fallback(pipe, X_tr, X_te)

def get_contribution_permutation_fallback(pipe, X_tr, X_te):
    """Fallback permutation method for SVM."""
    imp    = pipe.named_steps['imp']
    X_tr_imp = imp.transform(X_tr)
    train_med = np.nanmedian(X_tr_imp, axis=0)
    p_base = pipe.predict_proba(X_te)[0, 1]
    contribs = np.zeros(n_feat)
    for j in range(n_feat):
        X_perm    = X_te.copy()
        X_perm[0, j] = train_med[j]
        p_perm    = pipe.predict_proba(X_perm)[0, 1]
        contribs[j] = p_base - p_perm
    return contribs

print("Computing per-patient LOOCV contributions ...")
contrib_mat = np.zeros((n_pat, n_feat))   # rows=patients, cols=features
prob_loocv  = np.zeros(n_pat)

loo = LeaveOneOut()
for fold, (tr_idx, te_idx) in enumerate(loo.split(X_raw)):
    X_tr, X_te = X_raw[tr_idx], X_raw[te_idx]
    y_tr       = y[tr_idx]
    i          = te_idx[0]

    if is_ensemble:
        # Ensemble: average contribution across top3 models
        c_sum = np.zeros(n_feat)
        p_sum = 0.0
        for nm in top3_names:
            mdl = build_model(nm)
            if mdl is None: continue
            if 'GB' in nm:
                sw_tr = compute_sample_weight('balanced', y_tr)
                mdl.fit(X_tr, y_tr, clf__sample_weight=sw_tr)
            else:
                mdl.fit(X_tr, y_tr)
            p_sum += mdl.predict_proba(X_te)[0, 1]
            c_sum += get_contribution_shap(mdl, X_tr, X_te, nm)
        contrib_mat[i] = c_sum / len(top3_names)
        prob_loocv[i]  = p_sum / len(top3_names)
    else:
        nm  = model_type
        mdl = build_model(nm)
        if nm in ('GB d=2','LGBM'):
            sw_tr = compute_sample_weight('balanced', y_tr)
            mdl.fit(X_tr, y_tr, clf__sample_weight=sw_tr)
        else:
            mdl.fit(X_tr, y_tr)
        prob_loocv[i] = mdl.predict_proba(X_te)[0, 1]
        contrib_mat[i] = get_contribution_shap(mdl, X_tr, X_te, nm)

    if (fold + 1) % 10 == 0:
        print(f"  {fold+1}/{n_pat} done")

print(f"✓ Done. P(PR) corr with best_yp: {np.corrcoef(prob_loocv, best_yp)[0,1]:.3f}")

# ── 4. สร้าง tracker_contrib DataFrame ───────────────────────────────────
contrib_df = pd.DataFrame(contrib_mat, columns=feat_labels)
contrib_df.insert(0, 'pid',  pids)
contrib_df.insert(1, 'CT',   cts)
contrib_df.insert(2, 'prob', best_yp.round(3))
contrib_df.insert(3, 'pred', ['PR' if p > 0.5 else 'non-PR' for p in best_yp])
contrib_df.insert(4, 'ok',   [p == ('PR' if ct == 'PR' else 'non-PR')
                               for p, ct in zip(contrib_df['pred'], cts)])

# sort เหมือน tracker: PR → SD → PD
ct_ord = {'PR': 0, 'SD': 1, 'PD': 2}
contrib_df['_o'] = contrib_df['CT'].map(ct_ord)
contrib_df = contrib_df.sort_values(['_o', 'pid']).drop(columns='_o').reset_index(drop=True)

print(f"\ncontrib_df shape: {contrib_df.shape}")

# ── 5A. FIGURE 1: Summary Heatmap — Patient × Feature ─────────────────────
C_mat = contrib_df[feat_labels].values
vmax  = np.percentile(np.abs(C_mat), 95)
norm  = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

fig_h = max(8, n_pat * 0.35)
fig, ax = plt.subplots(figsize=(max(14, n_feat * 0.9), fig_h))
fig.patch.set_facecolor('white')

im = ax.imshow(C_mat, aspect='auto', cmap='RdBu_r', norm=norm)

# axis labels
ax.set_xticks(range(n_feat))
ax.set_xticklabels(feat_labels, rotation=45, ha='right', fontsize=8.5)
ax.set_yticks(range(len(contrib_df)))

# y-tick: PID + CT + prob + pred + ✓/✗
ylabels = []
for _, row in contrib_df.iterrows():
    icon = '✓' if row['ok'] else '✗'
    ylabels.append(f"{row['pid']}  [{row['CT']}]  P={row['prob']:.2f} {icon}")
ax.set_yticklabels(ylabels, fontsize=8.5)

# CT group separator lines
prev_ct = None
for ri, ct in enumerate(contrib_df['CT']):
    if ct != prev_ct and ri > 0:
        ax.axhline(ri - 0.5, color='black', lw=1.8)
    prev_ct = ct

# Colorbar
cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.02)
cb.set_label('Feature contribution\n(↑ pushes P(PR) up  |  ↓ pushes P(PR) down)',
             fontsize=9)

ax.set_title(
    f'Per-Patient Feature Contributions — {best_lbl}  (AUC={best_auc_f:.3f})\n'
    f'Red = pushes P(PR) ↑ (toward PR)   |   Blue = pushes P(PR) ↓ (toward non-PR)',
    fontsize=11, fontweight='bold', pad=12)

add_panel_labels(axes)
plt.tight_layout(pad=1.8)
plt.savefig(f'{OUT}figH1_contrib_heatmap.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.show()
print("✓ figH1_contrib_heatmap.png")

# ── 5B. FIGURE 2: Global Feature Importance (mean |contribution|) ──────────
mean_abs = np.abs(C_mat).mean(axis=0)
mean_pr  = C_mat[contrib_df['CT'] == 'PR'].mean(axis=0)
mean_npR = C_mat[contrib_df['CT'] != 'PR'].mean(axis=0)

order = np.argsort(mean_abs)[::-1]
feat_ord  = [feat_labels[i] for i in order]
ma_ord    = mean_abs[order]
mPR_ord   = mean_pr[order]
mNPR_ord  = mean_npR[order]

fig, axes = plt.subplots(1, 2, figsize=(16, max(5, n_feat * 0.35)))
fig.patch.set_facecolor('white')
fig.suptitle(f'Global Feature Importance and Directional Contribution — {best_lbl}',
             fontsize=12, fontweight='bold')

# Left: mean |contribution| bar
ax = axes[0]
bars = ax.barh(range(n_feat), ma_ord[::-1], color='#4393C3', alpha=0.85,
               edgecolor='white', linewidth=0.5)
ax.set_yticks(range(n_feat))
ax.set_yticklabels(feat_ord[::-1], fontsize=9)
ax.set_xlabel('Mean |contribution|', fontsize=10)
ax.set_title('Overall Feature Impact (Mean |Contribution|)', fontsize=11, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', alpha=0.25)

# Right: mean contribution PR vs non-PR (direction)
ax2 = axes[1]
y_pos = np.arange(n_feat)
ax2.barh(y_pos - 0.2, mPR_ord[::-1],  height=0.38, color='#2166AC',
         alpha=0.85, label='PR patients', edgecolor='white')
ax2.barh(y_pos + 0.2, mNPR_ord[::-1], height=0.38, color='#D6604D',
         alpha=0.85, label='non-PR patients', edgecolor='white')
ax2.axvline(0, color='black', lw=1.0)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(feat_ord[::-1], fontsize=9)
ax2.set_xlabel('Mean contribution\n(+) toward PR   (−) toward non-PR', fontsize=10)
ax2.set_title('Directional Contribution: PR versus Non-PR', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9, loc='lower right')
ax2.spines[['top', 'right']].set_visible(False)
ax2.grid(axis='x', alpha=0.25)

plt.tight_layout(pad=2.5)
plt.savefig(f'{OUT}figH2_global_importance.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.show()
print("✓ figH2_global_importance.png")

# ── 5C. FIGURE 3: Waterfall per patient ───────────────────────────────────
N_COLS   = 4
n_rows_w = int(np.ceil(n_pat / N_COLS))
TOP_N    = min(n_feat, 10)   # แสดง top-10 features ต่อคน

fig3, axes3 = plt.subplots(n_rows_w, N_COLS,
                            figsize=(N_COLS * 5.5, n_rows_w * 4.5))
fig3.patch.set_facecolor('white')
fig3.suptitle(
    f'Per-Patient SHAP Waterfall — {best_lbl}  (AUC = {best_auc_f:.3f})\n'
    f'Red bars = push P(PR) ↑   |   Blue bars = push P(PR) ↓   |   ✓ correct  ✗ wrong',
    fontsize=12, fontweight='bold', y=1.01)

for pi, (_, row) in enumerate(contrib_df.iterrows()):
    r  = pi // N_COLS
    c  = pi  % N_COLS
    ax = axes3[r, c] if n_rows_w > 1 else axes3[c]

    contribs = np.array([row[fl] for fl in feat_labels])
    order_p  = np.argsort(np.abs(contribs))[::-1][:TOP_N]
    c_sel    = contribs[order_p]
    f_sel    = [feat_labels[j] for j in order_p]

    colors_b = ['#D6604D' if v > 0 else '#4393C3' for v in c_sel]
    y_pos    = np.arange(len(f_sel))
    ax.barh(y_pos, c_sel[::-1], color=colors_b[::-1], alpha=0.88,
            edgecolor='white', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(f_sel[::-1], fontsize=7.5)
    ax.axvline(0, color='black', lw=0.8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    ax.tick_params(axis='x', labelsize=7)

    icon  = '✓' if row['ok'] else '✗'
    color = '#1B7837' if row['ok'] else '#762A83'
    ax.set_title(f"{row['pid']}  [{row['CT']}]  P={row['prob']:.2f}  {icon}",
                 fontsize=9, fontweight='bold', color=color, pad=4)


for pi in range(n_pat, n_rows_w * N_COLS):
    r = pi // N_COLS; c = pi % N_COLS
    ax_e = axes3[r, c] if n_rows_w > 1 else axes3[c]
    ax_e.set_visible(False)

plt.tight_layout(pad=2.0)
plt.savefig(f'{OUT}figH3_waterfall_all.png', dpi=160, bbox_inches='tight', facecolor='white')
plt.show()
print("✓ figH3_waterfall_all.png")

# ── 5D. Numeric summary table ─────────────────────────────────────────────
print("\n" + "="*70)
print(f"  FEATURE CONTRIBUTION SUMMARY — {best_lbl}")
print("="*70)
print(f"  {'Feature':<30}  {'Mean|C|':>8}  {'Mean PR':>9}  {'Mean non-PR':>11}")
print("  " + "-"*62)
for j in order:
    print(f"  {feat_labels[j]:<30}  {mean_abs[j]:>8.4f}  "
          f"{mean_pr[j]:>9.4f}  {mean_npR[j]:>11.4f}")
print("="*70)

# save contrib_df for downstream use
print(f"\ncontrib_df saved in memory. Shape: {contrib_df.shape}")
display(contrib_df.style
    .background_gradient(subset=feat_labels, cmap='RdBu_r',
                         vmin=-vmax, vmax=vmax)
    .format({fl: '{:+.3f}' for fl in feat_labels} |
            {'prob': '{:.3f}'})
    .set_caption(f'Per-Patient Feature Contributions — {best_lbl}'))

# %% [cell 57]
# ### Step 5.8 · Progression-Free Survival — Kaplan-Meier Analysis
#
# **Objective:** Validate model predictions against actual clinical outcomes using PFS.
#
# Requires `clinicalData.xlsx` with columns: `PID`, `PFS last check 7/18/2025`, `Status`.
#
# **Panels:**
# - **A:** PFS by CT response
# - **B:** PFS by model prediction
# - **C:** Correct vs misclassified patients
# - **D:** P(PR) score groups (>=0.5 vs <0.5)

# %% [cell 58]
# ── Step 5.8 · Progression-Free Survival (Kaplan-Meier) ────────────────────
# เชื่อมผลโมเดล (tracker) กับข้อมูล PFS จาก clinicalData.xlsx
# แสดง 4 กราฟ: A) CT Response  B) Model Prediction  C) Correct vs Misclassified  D) P(PR) score

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from scipy.stats import chi2
import pandas as pd

# Clinical data was loaded in the data-loading section.

# ── Merge tracker กับ PFS ──────────────────────────────────────────────────
surv = tracker[['pid', 'CT', 'prob', 'pred', 'ok']].copy()
surv = surv.merge(
    clinical[['PID', 'PFS last check 7/18/2025', 'Status']],
    left_on='pid', right_on='PID', how='inner'
)
surv['PFS_months'] = surv['PFS last check 7/18/2025'] / 30.44
surv['event']      = (~surv['Status'].str.strip().eq('Censored')).astype(int)
surv['PPR_group']  = surv['prob'].apply(lambda x: 'P(PR) ≥ 0.5' if x >= 0.5 else 'P(PR) < 0.5')
surv['match_lbl']  = surv['ok'].map({True: 'Correct', False: 'Misclassified'})
surv = surv.dropna(subset=['PFS_months']).reset_index(drop=True)

print(f"Patients with PFS data: {len(surv)}")
print(surv.groupby('CT')[['event']].agg(['count','sum']).to_string())

# ── KM functions ───────────────────────────────────────────────────────────
def km_with_ci(times, events):
    order = np.argsort(times)
    t = np.array(times, float)[order]
    e = np.array(events, int)[order]
    surv_val, gw = 1.0, 0.0
    ts, ss, lo_l, hi_l = [0.0], [1.0], [1.0], [1.0]
    for t_u in np.unique(t):
        d = e[t == t_u].sum()
        n = (t >= t_u).sum()
        if n > 0 and d > 0:
            surv_val *= (1 - d / n)
            gw += d / (n * (n - d)) if n > d else 0
        se = surv_val * np.sqrt(gw)
        ts.append(t_u); ss.append(surv_val)
        lo_l.append(max(0, surv_val - 1.96*se))
        hi_l.append(min(1, surv_val + 1.96*se))
    return np.array(ts), np.array(ss), np.array(lo_l), np.array(hi_l)

def log_rank_p(t1, e1, t2, e2):
    ev_t = np.unique(np.concatenate([t1[np.array(e1,bool)], t2[np.array(e2,bool)]]))
    O1, E1, V = 0, 0, 0
    for t in ev_t:
        n1=(t1>=t).sum(); n2=(t2>=t).sum()
        d1=((t1==t)&np.array(e1,bool)).sum(); d2=((t2==t)&np.array(e2,bool)).sum()
        n=n1+n2; d=d1+d2
        if n < 2: continue
        E1 += n1*d/n; O1 += d1
        if d < n: V += n1*n2*d*(n-d)/(n**2*(n-1))
    if V == 0: return 1.0
    return 1 - chi2.cdf((O1-E1)**2/V, df=1)

def s_at_t(ts, ss, t):
    return ss[max(0, np.searchsorted(ts, t, side='right')-1)]

def median_s(ts, ss):
    for i,s in enumerate(ss):
        if s <= 0.5: return round(ts[i],1)
    return None

def fmt_p(p):
    return 'P < 0.001' if p < 0.001 else f'P = {p:.3f}'

# ── Panel draw function ────────────────────────────────────────────────────
def draw_km(ax, groups, colors, order, pval_pairs, title,
            ylabel='Progression-free survival (proportion)',
            xmax=40, xticks=None, pval_x=0.38, pval_y=0.56):

    gdata = {}
    for lbl in order:
        g = groups[lbl]
        ts, ss, lo, hi = km_with_ci(g['PFS_months'].values, g['event'].values)
        gdata[lbl] = (ts, ss, lo, hi, g)

    # CI shading
    for lbl in order:
        ts, ss, lo, hi, _ = gdata[lbl]
        ax.fill_between(ts, lo, hi, step='post', alpha=0.13, color=colors[lbl])

    # KM lines + censored ticks
    for lbl in order:
        ts, ss, lo, hi, g = gdata[lbl]
        n_tot = len(g); n_ev = g['event'].sum()
        med = median_s(ts, ss)
        med_str = f'{med} mo' if med else 'NR'
        ax.step(ts, ss, where='post', color=colors[lbl], lw=2.2, zorder=3,
                label=f'{lbl}  (n={n_tot}, ev={n_ev}, mPFS={med_str})')
        c_t = g[g['event']==0]['PFS_months'].values
        c_s = [s_at_t(ts, ss, ct) for ct in c_t]
        ax.scatter(c_t, c_s, marker='+', color=colors[lbl], s=70, lw=1.8, zorder=5)

    # Log-rank p-values
    p_lines = []
    for l1, l2 in pval_pairs:
        g1=groups[l1]; g2=groups[l2]
        p = log_rank_p(g1['PFS_months'].values, g1['event'].values,
                       g2['PFS_months'].values, g2['event'].values)
        p_lines.append(f'{l1} vs {l2}   {fmt_p(p)}')
    ax.text(pval_x, pval_y, '\n'.join(p_lines), transform=ax.transAxes,
            fontsize=9, va='top',
            bbox=dict(boxstyle='round,pad=0.45', fc='white', ec='#cccccc', alpha=0.95))

    # Median reference line
    ax.axhline(0.5, color='#bbbbbb', ls=':', lw=1.0, zorder=1)

    # Formatting
    ax.set_xlim(-0.5, xmax); ax.set_ylim(-0.04, 1.08)
    ax.set_xlabel('Time (months)', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left', pad=10)
    if xticks: ax.set_xticks(xticks)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,_: f'{int(y*100)}%'))
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(0.8)
    ax.tick_params(direction='out', length=4, labelsize=10)
    ax.grid(axis='y', alpha=0.22, lw=0.6, ls='--')
    ax.legend(fontsize=8.5, loc='upper right', frameon=True,
              framealpha=0.9, edgecolor='#cccccc',
              handlelength=2.2, labelspacing=0.55)

# ── Define groups & colors ─────────────────────────────────────────────────
CT_MAP    = {'PR': 'PR', 'SD': 'SD', 'PD': 'PD'}
CT_COLORS = {'PR': '#2166AC', 'SD': '#4DAC26', 'PD': '#D6604D'}
ct_groups = {CT_MAP[k]: surv[surv['CT']==k] for k in ['PR','SD','PD']}

PRED_COLORS = {'Pred PR': '#2166AC', 'Pred non-PR': '#D6604D'}
pred_groups = {
    'Pred PR':     surv[surv['pred']=='PR'],
    'Pred non-PR': surv[surv['pred']=='non-PR'],
}

CORR_COLORS = {'Correct': '#1B7837', 'Misclassified': '#762A83'}
corr_groups = {
    'Correct':       surv[surv['match_lbl']=='Correct'],
    'Misclassified': surv[surv['match_lbl']=='Misclassified'],
}

PPR_COLORS = {'P(PR) ≥ 0.5': '#E08214', 'P(PR) < 0.5': '#4393C3'}
ppr_groups = {
    'P(PR) ≥ 0.5': surv[surv['PPR_group']=='P(PR) ≥ 0.5'],
    'P(PR) < 0.5': surv[surv['PPR_group']=='P(PR) < 0.5'],
}

# ── Plot: 4 panels ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.patch.set_facecolor('white')
fig.suptitle(
    f'Progression-Free Survival Analysis — cfDNA NSCLC  |  Model: {best_lbl}  (AUC = {best_auc_f:.3f})',
    fontsize=13, fontweight='bold', y=1.01
)

KW = dict(xmax=40, xticks=[0,10,20,30,40])

draw_km(axes[0,0], ct_groups, CT_COLORS,
        order=['PR','SD','PD'],
        pval_pairs=[('PR','SD'),('SD','PD')],
        title='A  PFS by CT Response', **KW)

draw_km(axes[0,1], pred_groups, PRED_COLORS,
        order=['Pred PR','Pred non-PR'],
        pval_pairs=[('Pred PR','Pred non-PR')],
        title=f'B  PFS by Model Prediction ({best_lbl})', **KW)

draw_km(axes[1,0], corr_groups, CORR_COLORS,
        order=['Correct','Misclassified'],
        pval_pairs=[('Correct','Misclassified')],
        title='C  PFS: Correctly vs Misclassified Patients', **KW)

draw_km(axes[1,1], ppr_groups, PPR_COLORS,
        order=['P(PR) ≥ 0.5','P(PR) < 0.5'],
        pval_pairs=[('P(PR) ≥ 0.5','P(PR) < 0.5')],
        title='D  PFS by P(PR) Score  (cutoff = 0.5)', **KW)

add_panel_labels(axes)
plt.tight_layout(pad=3.0)
fig_path = f'{OUT}/figG_km_pfs.png'
plt.savefig(fig_path, dpi=180, bbox_inches='tight', facecolor='white')
plt.show()
print(f"\n✓ Saved → {fig_path}")

# %% [cell 59]
# ── Step 5.8 Extra · PFS of Top-5 Models from Step 5.3 ───────────────────

# ── 1. เลือก top-5 จาก sweep (เรียงตาม AUC) ──────────────────────────────
sweep_sorted = sorted(sweep, key=lambda r: r['AUC'], reverse=True)

# กรอง unique model+k (ไม่ซ้ำ) แล้วเอา top-5
seen = set()
top5 = []
for r in sweep_sorted:
    key = (r['model'], r['k'])
    if key not in seen:
        seen.add(key)
        top5.append(r)
    if len(top5) == 5:
        break

print('Top-5 models selected:')
for i, r in enumerate(top5, 1):
    print(f"  {i}. {r['model']:12}  k={r['k']:2d}  AUC={r['AUC']:.3f}")

# ── 2. Map probability ของแต่ละโมเดลเข้า surv ────────────────────────────
for r in top5:
    col = f"top5_{r['model'].replace(' ','_').replace('=','')}_k{r['k']}"
    r['surv_col'] = col
    prob_map = dict(zip(df2['pid'].values, r['yp']))
    surv[col] = surv['pid'].map(prob_map)

# ── 3. สร้าง Figure: 2 rows × 3 cols (5 panels + 1 summary) ──────────────
fig, axes = plt.subplots(2, 3, figsize=(20, 11))
fig.patch.set_facecolor('white')
fig.suptitle(
    'Figure 20  |  Progression-Free Survival: Top-5 Models from Feature Selection\n'
    '(LOOCV | PR vs SD+PD, Ranked by AUC)',
    fontsize=14, fontweight='bold', y=1.01
)

RANK_COLORS = ['#C0392B', '#2980B9', '#27AE60', '#8E44AD', '#E67E22']
KM_COLORS   = {'Model PR': '#2166AC', 'Model non-PR': '#D6604D'}

for idx, r in enumerate(top5):
    row, col_ax = divmod(idx, 3)
    ax = axes[row][col_ax]
    col = r['surv_col']

    sub = surv.dropna(subset=[col]).copy()
    grp_pr   = sub[sub[col] >= 0.5]
    grp_npr  = sub[sub[col] <  0.5]

    if len(grp_pr) == 0 or len(grp_npr) == 0:
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                transform=ax.transAxes, fontsize=11)
        ax.set_visible(True)
        continue

    p_val = log_rank_p(
        grp_pr['PFS_months'].values,  grp_pr['event'].values,
        grp_npr['PFS_months'].values, grp_npr['event'].values
    )

    rank_color = RANK_COLORS[idx]
    title = (f'#{idx+1}  {r["model"]}  k={r["k"]}\n'
             f'AUC={r["AUC"]:.3f}  Acc={r["Acc"]:.3f}  '
             f'Sens={r["Sens"]:.3f}  Spec={r["Spec"]:.3f}')

    draw_km(ax,
            groups={'Model PR': grp_pr, 'Model non-PR': grp_npr},
            colors=KM_COLORS,
            order=['Model PR', 'Model non-PR'],
            pval_pairs=[('Model PR', 'Model non-PR')],
            title=title,
            xmax=40, xticks=[0, 10, 20, 30, 40],
            pval_x=0.36, pval_y=0.54)

    # กรอบสีตาม rank
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(rank_color)
        spine.set_linewidth(2.2)

# ── Panel 6 (row=1, col=2): Summary bar chart เปรียบเทียบ p-value ────────
ax_sum = axes[1][2]
ax_sum.set_facecolor('#fafafa')

model_names, p_vals, auc_vals, colors_bar = [], [], [], []
for idx, r in enumerate(top5):
    col = r['surv_col']
    sub = surv.dropna(subset=[col]).copy()
    g_pr  = sub[sub[col] >= 0.5]
    g_npr = sub[sub[col] <  0.5]
    if len(g_pr) == 0 or len(g_npr) == 0:
        continue
    p = log_rank_p(g_pr['PFS_months'].values,  g_pr['event'].values,
                   g_npr['PFS_months'].values, g_npr['event'].values)
    lbl = f"#{idx+1} {r['model']} k={r['k']}"
    model_names.append(lbl)
    p_vals.append(p)
    auc_vals.append(r['AUC'])
    colors_bar.append(RANK_COLORS[idx])

# เพิ่ม CT reference
ct_p = log_rank_p(
    surv[surv['CT_binary']=='PR']['PFS_months'].values,
    surv[surv['CT_binary']=='PR']['event'].values,
    surv[surv['CT_binary']=='non-PR']['PFS_months'].values,
    surv[surv['CT_binary']=='non-PR']['event'].values
)
model_names.append('CT Response\n(reference)')
p_vals.append(ct_p)
auc_vals.append(None)
colors_bar.append('#888888')

import matplotlib.patches as mpatches

y_pos = np.arange(len(model_names))
neg_log_p = [-np.log10(max(p, 1e-4)) for p in p_vals]

bars = ax_sum.barh(y_pos, neg_log_p, color=colors_bar, alpha=0.85,
                   edgecolor='white', linewidth=0.8)

# เส้น significance threshold (p=0.05)
ax_sum.axvline(-np.log10(0.05), color='#e74c3c', lw=1.6,
               ls='--', label='p = 0.05', zorder=5)

# label p-value บน bar
for i, (bar, p) in enumerate(zip(bars, p_vals)):
    p_str = 'p<0.001' if p < 0.001 else f'p={p:.3f}'
    ax_sum.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                p_str, va='center', ha='left', fontsize=9, fontweight='bold',
                color='#333333')

ax_sum.set_yticks(y_pos)
ax_sum.set_yticklabels(model_names, fontsize=9)
ax_sum.set_xlabel('−log₁₀(p-value)  [higher = more significant]', fontsize=10)
ax_sum.set_title('Summary of PFS Group Separation\n(Log-rank p-value)',
                 fontsize=11, fontweight='bold', loc='left', pad=8)
ax_sum.legend(fontsize=9, loc='lower right')
ax_sum.spines[['top', 'right']].set_visible(False)
ax_sum.grid(axis='x', alpha=0.25, lw=0.6)
ax_sum.set_xlim(0, max(neg_log_p) * 1.35)

# ── AUC annotation ด้านขวา ────────────────────────────────────────────────
ax_auc = ax_sum.twiny()
ax_auc.set_xlim(ax_sum.get_xlim())
ax_auc.set_xticks([])
for i, auc in enumerate(auc_vals):
    if auc is not None:
        ax_auc.text(max(neg_log_p) * 1.28, i,
                    f'AUC={auc:.3f}', va='center', ha='right',
                    fontsize=8, color='#555555', style='italic')

add_panel_labels(axes)
plt.tight_layout(pad=2.5)
plt.savefig(f'{OUT}figG4_km_top5_step53.png', dpi=180,
            bbox_inches='tight', facecolor='white')
plt.show()
print('✓ figG4_km_top5_step53.png')

# ── 4. Summary table ───────────────────────────────────────────────────────
print('\n' + '='*72)
print('  TOP-5 STEP 5.3 · PFS SEPARATION SUMMARY')
print('='*72)
print(f"  {'Rank':<5}  {'Model':12}  {'k':>3}  {'AUC':>6}  "
      f"{'n PR':>5}  {'n non-PR':>9}  {'mPFS PR':>8}  {'mPFS nPR':>9}  {'p-value':>10}")
print('  ' + '-'*68)

for idx, r in enumerate(top5):
    col = r['surv_col']
    sub = surv.dropna(subset=[col]).copy()
    g_pr  = sub[sub[col] >= 0.5]
    g_npr = sub[sub[col] <  0.5]
    if len(g_pr) == 0 or len(g_npr) == 0: continue

    p = log_rank_p(g_pr['PFS_months'].values,  g_pr['event'].values,
                   g_npr['PFS_months'].values, g_npr['event'].values)

    ts_pr, ss_pr, _, _ = km_with_ci(g_pr['PFS_months'].values,  g_pr['event'].values)
    ts_np, ss_np, _, _ = km_with_ci(g_npr['PFS_months'].values, g_npr['event'].values)
    med_pr  = median_s(ts_pr, ss_pr)
    med_npr = median_s(ts_np, ss_np)
    med_pr_s  = f'{med_pr} mo'  if med_pr  else 'NR'
    med_npr_s = f'{med_npr} mo' if med_npr else 'NR'

    sig = ' ✓' if p < 0.05 else ''
    print(f"  #{idx+1:<4}  {r['model']:12}  {r['k']:3d}  {r['AUC']:6.3f}  "
          f"{len(g_pr):5d}  {len(g_npr):9d}  "
          f"{med_pr_s:>8}  {med_npr_s:>9}  {fmt_p(p):>10}{sig}")

# CT reference
print('  ' + '-'*68)
g_pr_ct  = surv[surv['CT_binary']=='PR']
g_npr_ct = surv[surv['CT_binary']=='non-PR']
ts_pr, ss_pr, _, _ = km_with_ci(g_pr_ct['PFS_months'].values,  g_pr_ct['event'].values)
ts_np, ss_np, _, _ = km_with_ci(g_npr_ct['PFS_months'].values, g_npr_ct['event'].values)
ct_med_pr  = median_s(ts_pr, ss_pr); ct_med_npr = median_s(ts_np, ss_np)
print(f"  {'ref':<5}  {'CT Response':12}  {'–':>3}  {'–':>6}  "
      f"{len(g_pr_ct):5d}  {len(g_npr_ct):9d}  "
      f"{(str(ct_med_pr)+' mo') if ct_med_pr else 'NR':>8}  "
      f"{(str(ct_med_npr)+' mo') if ct_med_npr else 'NR':>9}  {fmt_p(ct_p):>10}")
print('='*72)

# %% [cell 60]
# ---
# ## Output Summary
#
# This notebook saves generated tables and figures under `./outputs/`.
#
# Main outputs include:
# - EDA figures for cfDNA profile comparison
# - Fragment length and end-motif visualizations
# - Stage 1 NSCLC vs Healthy classification results
# - Stage 2 PR vs non-PR treatment response prediction results
# - SHAP-based model interpretation plots
# - Per-patient prediction tracker
# - Kaplan-Meier progression-free survival plots
