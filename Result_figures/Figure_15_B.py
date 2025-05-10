import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Load the CSV ---
df = pd.read_csv('/mnt/scratch/nazilak/Results/variant_comparison.csv')

# --- Filter for TP53 or RB1 and correct variant_class ---
target_genes = ['TP53', 'RB1']
target_classes = ['NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS']

filtered = df[
    df['gene'].isin(target_genes) &
    df['variant_class'].isin(target_classes)
]

# --- Identify unique variants by chr, pos, ref, alt only ---
variant_cols = ['chr', 'pos', 'ref', 'alt']

# Separate by workflow
gatk = filtered[filtered['workflow'] == 'GATK'].drop_duplicates(subset=variant_cols)
clc = filtered[filtered['workflow'] == 'CLC'].drop_duplicates(subset=variant_cols)

# Create sets for comparison
gatk_set = set(tuple(x) for x in gatk[variant_cols].values)
clc_set = set(tuple(x) for x in clc[variant_cols].values)

# Unique and shared variants (by chr, pos, ref, alt)
unique_gatk = gatk_set - clc_set
unique_clc = clc_set - gatk_set
shared = gatk_set & clc_set

# For each gene, count how many variants (by chr,pos,ref,alt) are in each group
def count_variants(df, variant_set, gene):
    return df[
        df['gene'] == gene
    ].apply(lambda row: tuple(row[variant_cols]) in variant_set, axis=1).sum()

genes = target_genes
counts = {
    'GATK unique variants': [count_variants(gatk, unique_gatk, g) for g in genes],
    'Matching variants': [count_variants(gatk, shared, g) for g in genes],
    'CLC unique variants': [count_variants(clc, unique_clc, g) for g in genes]
}

# --- Plotting ---
labels = genes
gatk_only = counts['GATK unique variants']
both = counts['Matching variants']
clc_only = counts['CLC unique variants']

x = np.arange(len(labels))
width = 0.4  # Slightly wider bars

fig, ax = plt.subplots(figsize=(6, 6))  # Narrower figure

# Colors: light blue, salmon, light green
color_gatk = '#8ecae6'   # Light blue
color_both = '#ffb4a2'   # Salmon
color_clc = '#90be6d'    # Light green

# Plot stacked bars WITHOUT black edges
p1 = ax.bar(x, gatk_only, width, label='GATK unique variants', color=color_gatk)
p2 = ax.bar(x, both, width, bottom=gatk_only, label='Matching variants', color=color_both)
p3 = ax.bar(x, clc_only, width, bottom=np.array(gatk_only)+np.array(both), label='CLC unique variants', color=color_clc)

# Add counts inside each box
for i in range(len(labels)):
    # GATK only
    if gatk_only[i] > 0:
        ax.text(x[i], gatk_only[i]/2, str(gatk_only[i]), ha='center', va='center', fontsize=12, color='black')
    # Both
    if both[i] > 0:
        ax.text(x[i], gatk_only[i] + both[i]/2, str(both[i]), ha='center', va='center', fontsize=12, color='black')
    # CLC only
    if clc_only[i] > 0:
        ax.text(x[i], gatk_only[i] + both[i] + clc_only[i]/2, str(clc_only[i]), ha='center', va='center', fontsize=12, color='black')

# Labels and title
ax.set_ylabel('Number of Variants', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.legend(fontsize=11)

# Clean up spines for a cleaner look
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Bring bars closer to x-axis
ax.margins(y=0.02)

plt.tight_layout()
plt.savefig("/mnt/scratch/nazilak/Results/TP53_RB1.png", dpi=120)
plt.show()

