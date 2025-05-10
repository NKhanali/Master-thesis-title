    
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
################################## FIGURE 16_A.1 ######################
# Load the filtered CSV
df = pd.read_csv('/mnt/scratch/nazilak/Results/47_48_49_50.csv')
df_gatk_pass = df[
    (df['workflow'] == 'GATK') &
    (df['sample_nr'].isin([47, 48])) &
    (df['filter_status'] == 'PASS')
]


# Create unique variant identifier: chr:pos:ref:alt
df_gatk_pass['variant_id'] = df_gatk_pass['chr'].astype(str) + ':' + df_gatk_pass['pos'].astype(str) + ':' + df_gatk_pass['ref'] + '>' + df_gatk_pass['alt']

# Split by sample
variants_47 = set(df_gatk_pass[df_gatk_pass['sample_nr'] == 47]['variant_id'])
variants_48 = set(df_gatk_pass[df_gatk_pass['sample_nr'] == 48]['variant_id'])


# Create Venn diagram
plt.figure(figsize=(6,6))
venn2([variants_47, variants_48], set_labels=('Sample 47', 'Sample 48'))
plt.title('GATK Variant Overlap Between Sample 47 and 48')
plt.tight_layout()

# Save to file (choose your path and format)
output_path = '/mnt/scratch/nazilak/Results/venn_47_48_pass.png'
plt.savefig(output_path, dpi=300)
print(f"Venn diagram saved to: {output_path}")

# Optionally show the plot
plt.show()




################################## FIGURE 16_A.2 ######################
# Load the filtered CSV
df = pd.read_csv('/mnt/scratch/nazilak/Results/47_48_49_50.csv')

# Filter for GATK and samples 47 & 48
df_gatk = df[(df['workflow'] == 'GATK') & (df['sample_nr'].isin([47, 48]))]



# Create unique variant identifier: chr:pos:ref:alt
df_gatk['variant_id'] = df_gatk['chr'].astype(str) + ':' + df_gatk['pos'].astype(str) + ':' + df_gatk['ref'] + '>' + df_gatk['alt']

# Split by sample
variants_47 = set(df_gatk[df_gatk['sample_nr'] == 47]['variant_id'])
variants_48 = set(df_gatk[df_gatk['sample_nr'] == 48]['variant_id'])


# Create Venn diagram
plt.figure(figsize=(6,6))
venn2([variants_47, variants_48], set_labels=('Sample 47', 'Sample 48'))
plt.title('GATK Variant Overlap Between Sample 47 and 48')
plt.tight_layout()

# Save to file (choose your path and format)
output_path = '/mnt/scratch/nazilak/Results/venn_47_48_all.png'
plt.savefig(output_path, dpi=300)
print(f"Venn diagram saved to: {output_path}")

# Optionally show the plot
plt.show()
