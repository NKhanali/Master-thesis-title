    
import csv
import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict
from matplotlib_venn import venn3
# Input and output file names
input_file = '/mnt/scratch/nazilak/Results/variant_comparison.csv'
output_file = '/mnt/scratch/nazilak/Results/61_62_63.csv'

# Sample numbers and workflows to filter
sample_numbers = {'61', '62', '63'}
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


####################### FIGURE 19_A1 and 19_B1 ###################


# Input and output paths
input_file = '/mnt/scratch/nazilak/Results/61_62_63.csv'
output_venn_gatk = '/mnt/scratch/nazilak/Results/61_62_63_gatk_all.png'
output_venn_clc = '/mnt/scratch/nazilak/Results/61_62_63_clc_all.png'

# Dictionary to hold variant sets per sample for each workflow
gatk_variants = defaultdict(set)
clc_variants = defaultdict(set)

# Read the CSV and group variants
with open(input_file, mode='r', newline='') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        sample = row['sample_nr']
        workflow = row['workflow']
        
        # Create a unique variant identifier from the columns 'chr', 'pos', 'ref', 'alt'
        variant_id = f"{row['chr']}_{row['pos']}_{row['ref']}_{row['alt']}"

        # Add to the appropriate workflow and sample set
        if workflow == 'GATK':
            gatk_variants[sample].add(variant_id)
        elif workflow == 'CLC':
            clc_variants[sample].add(variant_id)

# Function to plot and save Venn diagrams
def save_venn(variant_dict, title, output_file):
    # Ensure that all three samples are included in the Venn diagram, even if empty
    set61 = variant_dict.get('61', set())
    set62 = variant_dict.get('62', set())
    set63 = variant_dict.get('63', set())

    # Make sure the sets are not empty for the Venn diagram
    venn3([set61, set62, set63], ('Sample 61', 'Sample 62', 'Sample 63'))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"{title} saved to {output_file}")

# Save both diagrams
save_venn(gatk_variants, 'GATK Workflow Variant Overlap (ALL VARIANTS)', output_venn_gatk)
save_venn(clc_variants, 'CLC Workflow Variant Overlap (ALL VARIANTS)', output_venn_clc)

####################### FIGURE 19_A2 and 19_B2 ###################

# Additional: Filtered GATK (PASS) and CLC (annotated) Venn diagrams

# Paths for filtered Venn outputs
output_venn_gatk_pass = '/mnt/scratch/nazilak/Results/61_62_63_gatk_pass.png'
output_venn_clc_annotated = '/mnt/scratch/nazilak/Results/61_62_63_clc_annotated.png'

# Dictionaries to hold filtered sets
gatk_pass_variants = defaultdict(set)
clc_annotated_variants = defaultdict(set)

# Read the CSV again and filter specifically for GATK "PASS" and CLC "annotated"
with open(input_file, mode='r', newline='') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        sample = row['sample_nr']
        workflow = row['workflow']
        variant_id = f"{row['chr']}_{row['pos']}_{row['ref']}_{row['alt']}"

        if workflow == 'GATK' and row.get('filter_status', '').upper() == 'PASS':
            gatk_pass_variants[sample].add(variant_id)
        elif workflow == 'CLC' and row.get('filter_status', '').lower() == 'annotated':
            clc_annotated_variants[sample].add(variant_id)

# Save the additional filtered diagrams
save_venn(gatk_pass_variants, 'GATK Workflow Variant Overlap (PASS Only)', output_venn_gatk_pass)
save_venn(clc_annotated_variants, 'CLC Workflow Variant Overlap (Annotated Only)', output_venn_clc_annotated)