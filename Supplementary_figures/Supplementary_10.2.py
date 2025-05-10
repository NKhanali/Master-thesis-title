import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
from matplotlib.patches import Patch  # Required for custom legend
from scipy.stats import pearsonr, spearmanr

# Load the CSV file containing the variant data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")

# Filter the data to focus on GATK (ALL) and CLC (ALL) variants
gatk = df[df['workflow'] == 'GATK'].copy()
clc = df[(df['workflow'] == 'CLC') & (df['filter_status'] == 'annotated')].copy()

# Remove duplicates (if any)
gatk = gatk.drop_duplicates(subset=['chr', 'pos', 'ref', 'alt', 'sample_nr'])
clc = clc.drop_duplicates(subset=['chr', 'pos', 'ref', 'alt', 'sample_nr'])

samples = range(1, 64)
gatk_unique_counts, clc_unique_counts, match_counts = [], [], []
all_matches = pd.DataFrame()

for sample in samples:
    gatk_sample = gatk[gatk['sample_nr'] == sample]
    clc_sample = clc[clc['sample_nr'] == sample]
    matches = pd.merge(gatk_sample, clc_sample, on=['chr', 'pos', 'ref', 'alt'], how='inner')
    match_count = len(matches)
    match_counts.append(match_count)
    unique_gatk = len(gatk_sample) - match_count
    unique_clc = len(clc_sample) - match_count
    gatk_unique_counts.append(unique_gatk)
    clc_unique_counts.append(unique_clc)
    all_matches = pd.concat([all_matches, matches])



#################################### FIGURE 12.B. ###############################
# Venn diagram data
total_gatk_variants = len(gatk)
total_clc_variants = len(clc)
total_matching_variants = len(all_matches)
unique_gatk_variants = total_gatk_variants - total_matching_variants
unique_clc_variants = total_clc_variants - total_matching_variants

# Plotting
fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[5, 4])  

# Bar plot
ax1 = fig.add_subplot(gs[0, 0])
indices = np.arange(len(samples))
bar_width = 0.5

ax1.bar(indices, gatk_unique_counts, width=bar_width, color='lightblue', label='Unique GATK Variants')
ax1.bar(indices, clc_unique_counts, width=bar_width, bottom=gatk_unique_counts, color='lightgreen', label='Unique CLC Variants')
ax1.bar(indices, match_counts, width=bar_width, bottom=np.array(gatk_unique_counts)+np.array(clc_unique_counts), color='salmon', label='Matching Variants')

total_counts = np.array(gatk_unique_counts) + np.array(clc_unique_counts) + np.array(match_counts)
max_count = total_counts.max()
ax1.set_ylim(0, max_count + 10)

ax1.set_xticks(indices[::3])
ax1.set_xticklabels([str(samples[i]) for i in range(0, len(samples), 3)], rotation=90)
ax1.set_xlabel("Sample Number", fontsize=16)
ax1.set_ylabel("Number of Variants", fontsize=16)

# Remove the legend from ax1
# Custom legend handles for the bar plot
legend_handles = [
    Patch(color='lightblue', label='Unique GATK Variants'),
    Patch(color='lightgreen', label='Unique CLC Variants'),
    Patch(color='salmon', label='Matching Variants')
]

# Venn diagram
ax2 = fig.add_subplot(gs[0, 1])
venn = venn2(
    subsets=(unique_gatk_variants, unique_clc_variants, total_matching_variants),
    set_labels=('', ''),
    set_colors=('lightblue', 'lightgreen'),
    ax=ax2
)

venn.get_patch_by_id('11').set_color('salmon')
venn.get_patch_by_id('10').set_color('lightblue')
venn.get_patch_by_id('01').set_color('lightgreen')

# Move the legend above the Venn diagram
ax2.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=1, fontsize=12, frameon=False)

# Adjust layout and spacing
plt.tight_layout()
plt.subplots_adjust(wspace=0.3)

# Save the figure
plt.savefig("/mnt/scratch/nazilak/Results/Matching_variants/combined_matching_variants3.png")
plt.show()


########################### SUPPLEMENTARY 10_A2_B2_C2 ############################
#correlation between the number of unique variants (from both workflows) and the number of matches (overlaps) across your samples?


# === 1. Load Data ===
summary_df = pd.read_csv("/mnt/scratch/nazilak/Results/Matching_variants/matching_summary3.csv")

# === 2. Calculate sum of unique variants ===
summary_df['unique_sum'] = summary_df['unique_gatk'] + summary_df['unique_clc']

# === 3. Define plotting function ===
def scatter_with_regression(x, y, xlabel, ylabel, out_file):
    plt.figure(figsize=(8,6))
    plt.scatter(x, y, color='blue', alpha=0.7, label='Samples')
    # Regression line
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m*x + b, color='red', linestyle='--', label='Fit')
    # Correlation
    r, p = pearsonr(x, y)
    plt.title(f"{xlabel} vs {ylabel}\nPearson r={r:.2f}, p={p:.3g}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

# === 4. Make three scatter plots ===

# 1. Unique GATK + Unique CLC
scatter_with_regression(
    summary_df['unique_sum'],
    summary_df['match'],
    'Unique GATK + Unique CLC',
    'Matching Variants',
    "/mnt/scratch/nazilak/Results/Matching_variants/3scatter_sum_vs_match.png"
)

# 2. Unique GATK only
scatter_with_regression(
    summary_df['unique_gatk'],
    summary_df['match'],
    'Unique GATK',
    'Matching Variants',
    "/mnt/scratch/nazilak/Results/Matching_variants/3scatter_gatk_vs_match.png"
)

# 3. Unique CLC only
scatter_with_regression(
    summary_df['unique_clc'],
    summary_df['match'],
    'Unique CLC',
    'Matching Variants',
    "/mnt/scratch/nazilak/Results/Matching_variants/3scatter_clc_vs_match.png"
)

print("All scatter plots saved successfully.")
