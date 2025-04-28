import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr
# Load the CSV file
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")

# Filter the data
gatk_data = df[(df['workflow'] == 'GATK') & (df['filter_status'] == 'PASS')]
clc_data = df[(df['workflow'] == 'CLC') & (df['filter_status'] == 'annotated')]



# Count variants for each sample and workflow
gatk_counts = gatk_data.groupby('sample_nr').size().reset_index(name='count')
gatk_counts['workflow'] = 'GATK'
clc_counts = clc_data.groupby('sample_nr').size().reset_index(name='count')
clc_counts['workflow'] = 'CLC'

# Merge for correlation
merged_counts = pd.merge(
    gatk_counts[['sample_nr', 'count']],
    clc_counts[['sample_nr', 'count']],
    on='sample_nr',
    suffixes=('_gatk', '_clc')
)

print(gatk_counts)
print(clc_counts)
# Pearson correlation
pearson_corr, _ = pearsonr(merged_counts['count_gatk'], merged_counts['count_clc'])

# Total counts
total_gatk_variants = gatk_data.shape[0]
total_clc_variants = clc_data.shape[0]


# Combine counts
combined_counts = pd.concat([gatk_counts, clc_counts])

# Plot setup
plt.figure(figsize=(12, 6))

# Barplot in background
sns.barplot(x='sample_nr', y='count', hue='workflow', data=combined_counts, dodge=True, alpha=0.6)

# Stripplot overlay
sns.stripplot(x='sample_nr', y='count', hue='workflow', data=combined_counts, dodge=True, marker='o', edgecolor='auto', linewidth=0.5)

# Fix legend (due to double entries from strip + bar)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), title='Workflow')

# Customize
plt.title(
    f'Variants per Sample from GATK (PASS) and CLC (Annotated)\n'
    f'Pearson Correlation: {pearson_corr:.2f} | Total - GATK: {total_gatk_variants}, CLC: {total_clc_variants}',
    fontsize=16
)
plt.xlabel('Sample Number', fontsize=12)
plt.ylabel('Number of Variants', fontsize=12)
plt.xticks(rotation=90)
plt.tight_layout()

# Save
os.makedirs('/mnt/scratch/nazilak/Results/Number_of_variants/Default/', exist_ok=True)
plt.savefig('/mnt/scratch/nazilak/Results/Number_of_variants/Default/plot1_GATK(PASS)_CLC(ANN).png')

# Show
plt.show()
