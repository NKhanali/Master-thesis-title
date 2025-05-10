import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr

# Load the CSV file
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")

# Filter the data
gatk_data = df[df['workflow'] == 'GATK']# Include all GATK variants (no filter on 'filter_status')
clc_data = df[(df['workflow'] == 'CLC') & (df['filter_status'] == 'annotated')] # Include only annotated CLC variants 

# Count variants for each sample and workflow
gatk_counts = gatk_data.groupby('sample_nr').size().reset_index(name='count')
gatk_counts['workflow'] = 'GATK'
clc_counts = clc_data.groupby('sample_nr').size().reset_index(name='count')
clc_counts['workflow'] = 'CLC'

# Merge counts for comparison using outer join to include all samples
merged_counts = pd.merge(
    gatk_counts[['sample_nr', 'count']],
    clc_counts[['sample_nr', 'count']],
    on='sample_nr',
    how='outer',
    suffixes=('_gatk', '_clc')
).fillna(0)

# Ensure counts are integers
merged_counts['count_gatk'] = merged_counts['count_gatk'].astype(int)
merged_counts['count_clc'] = merged_counts['count_clc'].astype(int)

# Calculate difference and absolute difference
merged_counts['difference'] = merged_counts['count_gatk'] - merged_counts['count_clc']
merged_counts['abs_difference'] = merged_counts['difference'].abs()

# Combine counts for plotting
gatk_counts_renamed = gatk_counts.rename(columns={'count': 'count', 'sample_nr': 'sample_nr'})
clc_counts_renamed = clc_counts.rename(columns={'count': 'count', 'sample_nr': 'sample_nr'})
combined_counts = pd.concat([gatk_counts_renamed, clc_counts_renamed])

# Define custom palette
custom_palette = {'CLC': '#60c060', 'GATK': '#4a9cd6'}  # Darker shades

# Plotting
plt.figure(figsize=(7, 8))

# Barplot
sns.barplot(
    x='sample_nr',
    y='count',
    hue='workflow',
    data=combined_counts,
    dodge=True,
    alpha=0.7,
    palette=custom_palette
)

# Stripplot
sns.stripplot(
    x='sample_nr',
    y='count',
    hue='workflow',
    data=combined_counts,
    dodge=True,
    marker='o',
    edgecolor='auto',
    linewidth=0.5,
    palette=custom_palette
)

# Reduce x-axis clutter
plt.xticks(range(0, 63, 5))  # Show every 5th sample number

# Fix legend (remove duplicates)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(
    by_label.values(),
    by_label.keys(),
    title='Workflow',
    fontsize=14,
    title_fontsize=16,
    markerscale=2,
    loc='upper left',
    frameon=True,
    borderaxespad=1.0
)

# Averages
avg_gatk = gatk_counts['count'].mean()
avg_clc = clc_counts['count'].mean()

plt.axhline(avg_gatk, color='blue', linestyle='--', linewidth=2, label=f'GATK Avg: {avg_gatk:.0f}')
plt.axhline(avg_clc, color='green', linestyle='--', linewidth=2, label=f'CLC Avg: {avg_clc:.0f}')

# Update legend to include avg lines
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(
    by_label.values(),
    by_label.keys(),
    title='Workflow',
    fontsize=14,
    title_fontsize=16,
    markerscale=2,
    loc='upper left',
    frameon=True,
    borderaxespad=1.0
)

# Customize
plt.xlabel('Sample Number', fontsize=16)
plt.ylabel('Number of Variants', fontsize=16)
plt.xticks(rotation=90)
plt.tight_layout()

# Save
os.makedirs('/mnt/scratch/nazilak/Results/Number_of_variants/Default/', exist_ok=True)
plt.savefig('/mnt/scratch/nazilak/Results/Number_of_variants/Default/plot5_GATK(ALL)_CLC(ANN).png')

# Show
plt.show()
