import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
from matplotlib.patches import Patch  # Required for custom legend
from scipy.stats import pearsonr, spearmanr

# Load data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")
gatk = df[(df['workflow'] == 'GATK') & (df['filter_status'] == 'PASS')].copy()
clc = df[(df['workflow'] == 'CLC') & (df['filter_status'] == 'annotated')].copy()

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
gs = fig.add_gridspec(2, 2, height_ratios=[0.5, 5], width_ratios=[5, 3])

# Bar plot with broken y-axis
ax1_top = fig.add_subplot(gs[0, 0])
ax1_bottom = fig.add_subplot(gs[1, 0], sharex=ax1_top)

indices = np.arange(len(samples))
bar_width = 0.5

# Upper part (above 200)
ax1_top.bar(indices, gatk_unique_counts, width=bar_width, color='lightblue')
ax1_top.bar(indices, clc_unique_counts, width=bar_width, bottom=gatk_unique_counts, color='lightgreen')
ax1_top.bar(indices, match_counts, width=bar_width, bottom=np.array(gatk_unique_counts)+np.array(clc_unique_counts), color='salmon')
ax1_top.set_ylim(200, 270)
ax1_top.spines['bottom'].set_visible(False)
ax1_top.tick_params(labelbottom=False)

# Lower part (below 200)
ax1_bottom.bar(indices, gatk_unique_counts, width=bar_width, color='lightblue')
ax1_bottom.bar(indices, clc_unique_counts, width=bar_width, bottom=gatk_unique_counts, color='lightgreen')
ax1_bottom.bar(indices, match_counts, width=bar_width, bottom=np.array(gatk_unique_counts)+np.array(clc_unique_counts), color='salmon')
ax1_bottom.set_ylim(0, 200)
ax1_bottom.spines['top'].set_visible(False)

# Diagonal lines to indicate break
d = .005
kwargs = dict(transform=ax1_top.transAxes, color='k', clip_on=False)
ax1_top.plot((-d, +d), (-d, +d), **kwargs)
ax1_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

kwargs.update(transform=ax1_bottom.transAxes)
ax1_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax1_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

ax1_bottom.set_xticks(indices[::3])
ax1_bottom.set_xticklabels([str(samples[i]) for i in range(0, len(samples), 3)], rotation=90)
ax1_bottom.set_xlabel("Sample Number", fontsize=16)
ax1_bottom.set_ylabel("Number of Variants", fontsize=16)

# Venn diagram
ax2 = fig.add_subplot(gs[:, 1])
venn = venn2(
    subsets=(unique_gatk_variants, unique_clc_variants, total_matching_variants),
    set_labels=('', ''),
    set_colors=('lightblue', 'lightgreen'),
    ax=ax2
)
venn.get_patch_by_id('11').set_color('salmon')
venn.get_patch_by_id('10').set_color('lightblue')
venn.get_patch_by_id('01').set_color('lightgreen')

# Custom legend above Venn diagram
legend_handles = [
    Patch(color='lightblue', label='Unique GATK Variants'),
    Patch(color='lightgreen', label='Unique CLC Variants'),
    Patch(color='salmon', label='Matching Variants')
]
ax2.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=1, fontsize=12, frameon=False)

plt.tight_layout()
plt.subplots_adjust(hspace=0.05, wspace=0.3)
plt.savefig("/mnt/scratch/nazilak/Results/Matching_variants/combined_matching_variants.png")
plt.show()

