
import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2


input_file = '/mnt/scratch/nazilak/Results/variant_comparison.csv'
output_file = '/mnt/scratch/nazilak/Results/55_56_filtered.csv'

sample_numbers = {'55', '56'}
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
    print("No matching variants found.")
    exit()


df = pd.read_csv(output_file)

df['variant_id'] = df['chr'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + '>' + df['alt']


def get_variants(dataframe, sample, workflow, passed_only=False, annotated_only=False):
    subset = dataframe[
        (dataframe['sample_nr'] == int(sample)) &
        (dataframe['workflow'] == workflow)
    ]
    if passed_only:
        subset = subset[subset['filter_status'].str.upper() == 'PASS']
    if annotated_only:
        subset = subset[subset['filter_status'].str.lower() == 'annotated']
    return set(subset['variant_id'])



def plot_venn(set1, set2, labels, title, output_path):
    plt.figure(figsize=(6, 6))
    venn2([set1, set2], set_labels=labels)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")
    plt.close()


clc_55_pass_ann = get_variants(df, '55', 'CLC', passed_only=False, annotated_only=True)
clc_56_pass_ann = get_variants(df, '56', 'CLC', passed_only=False, annotated_only=True)
gatk_56_pass = get_variants(df, '56', 'GATK', passed_only=True)

plot_venn(clc_55_pass_ann, clc_56_pass_ann,
          ('CLC (55)', 'CLC (56)'), 'Passed/Annotated: CLC 55 vs 56',
          '/mnt/scratch/nazilak/Results/venn1_passed_clc_55_vs_56.png')

plot_venn(clc_55_pass_ann, gatk_56_pass,
          ('CLC (55)', 'GATK (56)'), 'Passed/Annotated: CLC 55 vs GATK 56',
          '/mnt/scratch/nazilak/Results/venn2_passed_clc_55_vs_gatk_56.png')

plot_venn(clc_56_pass_ann, gatk_56_pass,
          ('CLC (56)', 'GATK (56)'), 'Passed/Annotated: CLC 56 vs GATK 56',
          '/mnt/scratch/nazilak/Results/venn3_passed_clc_56_vs_gatk_56.png')


clc_55_all = get_variants(df, '55', 'CLC')
clc_56_all = get_variants(df, '56', 'CLC')
gatk_56_all = get_variants(df, '56', 'GATK')

plot_venn(clc_55_all, clc_56_all,
          ('CLC (55)', 'CLC (56)'), 'All Variants: CLC 55 vs 56',
          '/mnt/scratch/nazilak/Results/venn4_all_clc_55_vs_56.png')

plot_venn(clc_55_all, gatk_56_all,
          ('CLC (55)', 'GATK (56)'), 'All Variants: CLC 55 vs GATK 56',
          '/mnt/scratch/nazilak/Results/venn5_all_clc_55_vs_gatk_56.png')

plot_venn(clc_56_all, gatk_56_all,
          ('CLC (56)', 'GATK (56)'), 'All Variants: CLC 56 vs GATK 56',
          '/mnt/scratch/nazilak/Results/venn6_all_clc_56_vs_gatk_56.png')
  
