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


