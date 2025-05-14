    
import csv
import pandas as pd
from matplotlib import pyplot as plt
from collections import defaultdict
from matplotlib_venn import venn3

input_file = '/mnt/scratch/nazilak/Results/variant_comparison.csv'
output_file = '/mnt/scratch/nazilak/Results/61_62_63.csv'


sample_numbers = {'61', '62', '63'}
workflows = {'GATK', 'CLC'}

with open(input_file, mode='r', newline='') as infile:
    reader = csv.DictReader(infile)
    filtered_rows = [
        row for row in reader
        if row['sample_nr'] in sample_numbers and row['workflow'] in workflows
    ]

if filtered_rows:
    with open(output_file, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(filtered_rows)
    print(f"Filtered variants written to {output_file}")
else:
    print("No matching variants found for the specified filters.")


####################### FIGURE 19_A1 and 19_B1 ###################
input_file = '/mnt/scratch/nazilak/Results/61_62_63.csv'
output_venn_gatk = '/mnt/scratch/nazilak/Results/61_62_63_gatk_all.png'
output_venn_clc = '/mnt/scratch/nazilak/Results/61_62_63_clc_all.png'

gatk_variants = defaultdict(set)
clc_variants = defaultdict(set)

with open(input_file, mode='r', newline='') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        sample = row['sample_nr']
        workflow = row['workflow']
        

        variant_id = f"{row['chr']}_{row['pos']}_{row['ref']}_{row['alt']}"

        if workflow == 'GATK':
            gatk_variants[sample].add(variant_id)
        elif workflow == 'CLC':
            clc_variants[sample].add(variant_id)


def save_venn(variant_dict, title, output_file):
    set61 = variant_dict.get('61', set())
    set62 = variant_dict.get('62', set())
    set63 = variant_dict.get('63', set())


    venn3([set61, set62, set63], ('Sample 61', 'Sample 62', 'Sample 63'))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"{title} saved to {output_file}")

save_venn(gatk_variants, 'GATK Workflow Variant Overlap (ALL VARIANTS)', output_venn_gatk)
save_venn(clc_variants, 'CLC Workflow Variant Overlap (ALL VARIANTS)', output_venn_clc)

####################### FIGURE 19_A2 and 19_B2 ###################

output_venn_gatk_pass = '/mnt/scratch/nazilak/Results/61_62_63_gatk_pass.png'
output_venn_clc_annotated = '/mnt/scratch/nazilak/Results/61_62_63_clc_annotated.png'

gatk_pass_variants = defaultdict(set)
clc_annotated_variants = defaultdict(set)


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

save_venn(gatk_pass_variants, 'GATK Workflow Variant Overlap (PASS Only)', output_venn_gatk_pass)
save_venn(clc_annotated_variants, 'CLC Workflow Variant Overlap (Annotated Only)', output_venn_clc_annotated)