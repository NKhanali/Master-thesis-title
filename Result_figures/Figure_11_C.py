import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr


df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")


gatk_data = df[df['workflow'] == 'GATK']  
clc_data = df[df['workflow'] == 'CLC']    

gatk_counts = gatk_data.groupby('sample_nr').size().reset_index(name='count')
gatk_counts['workflow'] = 'GATK'
clc_counts = clc_data.groupby('sample_nr').size().reset_index(name='count')
clc_counts['workflow'] = 'CLC'

merged_counts = pd.merge(
    gatk_counts[['sample_nr', 'count']],
    clc_counts[['sample_nr', 'count']],
    on='sample_nr',
    how='outer',
    suffixes=('_gatk', '_clc')
).fillna(0)

merged_counts['count_gatk'] = merged_counts['count_gatk'].astype(int)
merged_counts['count_clc'] = merged_counts['count_clc'].astype(int)

merged_counts['difference'] = merged_counts['count_gatk'] - merged_counts['count_clc']
merged_counts['abs_difference'] = merged_counts['difference'].abs()

gatk_counts_renamed = gatk_counts.rename(columns={'count': 'count', 'sample_nr': 'sample_nr'})
clc_counts_renamed = clc_counts.rename(columns={'count': 'count', 'sample_nr': 'sample_nr'})
combined_counts = pd.concat([gatk_counts_renamed, clc_counts_renamed])

custom_palette = {'CLC': '#60c060', 'GATK': '#4a9cd6'}  


plt.figure(figsize=(7, 8))


sns.barplot(
    x='sample_nr',
    y='count',
    hue='workflow',
    data=combined_counts,
    dodge=True,
    alpha=0.7,
    palette=custom_palette
)

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


plt.xticks(rotation=90)
plt.xticks(range(0, 63, 5))  


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

avg_gatk = gatk_counts['count'].mean()
avg_clc = clc_counts['count'].mean()

plt.axhline(avg_gatk, color='blue', linestyle='--', linewidth=2, label=f'GATK Avg: {avg_gatk:.0f}')
plt.axhline(avg_clc, color='green', linestyle='--', linewidth=2, label=f'CLC Avg: {avg_clc:.0f}')

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

plt.xlabel('Sample Number', fontsize=16)
plt.ylabel('Number of Variants', fontsize=16)
plt.xticks(rotation=90)
plt.tight_layout()

os.makedirs('/mnt/scratch/nazilak/Results/Number_of_variants/Default/', exist_ok=True)
plt.savefig('/mnt/scratch/nazilak/Results/Number_of_variants/Default/plot3_GATK(ALL)_CLC(ALL).png')

plt.show()
