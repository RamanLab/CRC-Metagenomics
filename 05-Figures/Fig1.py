
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.ticker import LogLocator, LogFormatterSciNotation


CSV_PATH   = "FinalReadSummary.csv"   # path to your CSV file
SAMPLE_COL = "Sample_Name"        # column with sample identifiers
LABEL_COL  = "Labels"             # column with group labels (e.g. NT / T)

STEP_LABELS = [
    "Raw Reads",
    "Trimmed Reads",
    "Host-aligned Reads",
    "Host-unaligned Reads",
    "Duplicate removed host Reads",
    "Decontaminated host Reads",
    "Kraken unclassified Reads",
    "Bracken classified Reads",
    "Human Reads (Bracken)",
    "Microbial Reads (Bracken)"
]

COLORS = ["#a8c8e8", "#1a3a6b"]   # light blue -> dark navy

OUTPUT_PATH = "read_summary_boxplot_total_analysis.pdf"

# LOAD & VALIDATE DATA

df = pd.read_csv(CSV_PATH)

# Identify numeric (pipeline-step) columns automatically
exclude_cols = {SAMPLE_COL, LABEL_COL}
step_cols = [c for c in df.columns if c not in exclude_cols and
             pd.api.types.is_numeric_dtype(df[c])]

if not step_cols:
    raise ValueError("No numeric columns found. Check SAMPLE_COL and LABEL_COL settings.")

print("Length of steps cols:", len(step_cols))
print("Length of step labels:", len(STEP_LABELS))
# Use user-supplied labels or fall back to column names
if STEP_LABELS and len(STEP_LABELS) == len(step_cols):
    x_labels = STEP_LABELS
else:
    print("Warning: STEP_LABELS length doesn't match numeric columns; using column names.")
    x_labels = step_cols

n_steps  = len(step_cols)
groups   = sorted(df[LABEL_COL].unique())          # e.g. ['NT', 'T']
n_groups = len(groups)

if n_groups > len(COLORS):
    cmap   = plt.cm.get_cmap("Blues", n_groups + 2)
    COLORS = [cmap(i / (n_groups + 1)) for i in range(1, n_groups + 1)]

color_map = {g: COLORS[i] for i, g in enumerate(groups)}

# LAYOUT GEOMETRY 

width   = 0.30
gap     = 0.06
spacing = 1.1

total_span = (n_groups - 1) * (width + gap)
offsets    = [i * (width + gap) - total_span / 2 for i in range(n_groups)]

positions = {
    g: [i * spacing + offsets[gi] for i in range(n_steps)]
    for gi, g in enumerate(groups)
}

# DRAW PLOT

fig, ax = plt.subplots(figsize=(max(10, n_steps * 1.8), 7))

rng = np.random.default_rng(42)

def draw_boxes(ax, group, pos_list, data_list, color):
    ax.boxplot(
        data_list,
        positions=pos_list,
        widths=width,
        patch_artist=True,
        notch=False,
        showfliers=False,
        medianprops=dict(color="#e8c040", linewidth=2),
        whiskerprops=dict(color=color, linewidth=1.2),
        capprops=dict(color=color, linewidth=1.2),
        boxprops=dict(facecolor=color, color=color, alpha=0.85),
        zorder=2,
    )
    for pos, vals in zip(pos_list, data_list):
        jx = rng.uniform(-width * 0.45, width * 0.45, size=len(vals))
        ax.scatter(pos + jx, vals, color="#aaaaaa", s=10, alpha=0.55, zorder=3)

for group in groups:
    group_data = [df.loc[df[LABEL_COL] == group, col].values for col in step_cols]
    draw_boxes(ax, group, positions[group], group_data, color_map[group])

ax.set_yscale("log")

# AXES & FORMATTING 

mid_positions = [np.mean([positions[g][i] for g in groups]) for i in range(n_steps)]
ax.set_xticks(mid_positions)
ax.set_xticklabels(x_labels, rotation=30, ha="right", fontsize=13)

ax.set_ylabel("Number of Reads", fontsize=13)
ax.set_xlim(
    min(positions[groups[0]]) - 0.55,
    max(positions[groups[-1]]) + 0.55,
)

ax.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10))
ax.tick_params(axis="y", labelsize=9)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)

# LEGEND

patches = [
    mpatches.Patch(facecolor=color_map[g], edgecolor=color_map[g], label=g)
    for g in groups
]
ax.legend(handles=patches, title=LABEL_COL, fontsize=12,
          title_fontsize=10, loc="upper right", framealpha=0.8)

#plt.title("Read counts at each pipeline step", fontsize=13, pad=12)
plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
print(f"Saved -> {OUTPUT_PATH}")
plt.show()