    
import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from upsetplot import UpSet, from_memberships

# Input and output file names
input_file = '/mnt/scratch/nazilak/Results/variant_comparison.csv'
output_file = '/mnt/scratch/nazilak/Results/47_48_49_50.csv'

# Sample numbers and workflows to filter
sample_numbers = {'47', '48', '49', '50'}
workflows = {'GATK', 'CLC'}

# Open the input file for reading
with open(input_file, mode='r', newline='') as infile:
    reader = csv.DictReader(infile)
    filtered_rows = [
        row for row in reader
        if row['sample_nr'] in sample_numbers and row['workflow'] in workflows
    ]

# Write the filtered data to the output file
if filtered_rows:
    with open(output_file, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(filtered_rows)
    print(f"Filtered variants written to {output_file}")
else:
    print("No matching variants found for the specified filters.")

####################### FIGURE 17_A ##########################
# Filter for CLC and samples 49 & 50
df_clc = df[
    (df['workflow'] == 'GATK') &
    (df['sample_nr'].isin([49, 50]))
]

# Create unique variant identifier: chr:pos:ref:alt
df_clc['variant_id'] = df_clc['chr'].astype(str) + ':' + df_clc['pos'].astype(str) + ':' + df_clc['ref'] + '>' + df_clc['alt']

# Split by sample
variants_49 = set(df_clc[df_clc['sample_nr'] == 49]['variant_id'])
variants_50 = set(df_clc[df_clc['sample_nr'] == 50]['variant_id'])

# Create Venn diagram
plt.figure(figsize=(6,6))
venn2([variants_49, variants_50], set_labels=('Sample 49', 'Sample 50'))
plt.title('GATK Variant Overlap Between Sample 49 and 50')
plt.tight_layout()

# Save to file
output_path = '/mnt/scratch/nazilak/Results/49_50gatk_match.png'
plt.savefig(output_path, dpi=300)
print(f"Venn diagram saved to: {output_path}")

# Optionally show the plot
plt.show()


####################### FIGURE 17_B ##########################
# Filter for CLC and samples 49 & 50
df_clc = df[
    (df['workflow'] == 'CLC') &
    (df['sample_nr'].isin([49, 50]))
]

# Create unique variant identifier: chr:pos:ref:alt
df_clc['variant_id'] = df_clc['chr'].astype(str) + ':' + df_clc['pos'].astype(str) + ':' + df_clc['ref'] + '>' + df_clc['alt']

# Split by sample
variants_49 = set(df_clc[df_clc['sample_nr'] == 49]['variant_id'])
variants_50 = set(df_clc[df_clc['sample_nr'] == 50]['variant_id'])

# Create Venn diagram
plt.figure(figsize=(6,6))
venn2([variants_49, variants_50], set_labels=('Sample 49', 'Sample 50'))
plt.title('CLC Variant Overlap Between Sample 49 and 50')
plt.tight_layout()

# Save to file
output_path = '/mnt/scratch/nazilak/Results/49_50clc_match.png'
plt.savefig(output_path, dpi=300)
print(f"Venn diagram saved to: {output_path}")

# Optionally show the plot
plt.show()
