import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns

# Define the function to parse VCF files
def vcf(file_path, workflow, sample_type=None):
    data = []
    samples = []

    with open(file_path, 'r', encoding='latin1') as f:
        for line in f:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = list(set([sample.split('_')[0] for sample in parts[9:]]))
            elif not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom = parts[0]
                if workflow == 'CLC' and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                info = parts[7]

                if workflow == 'GATK':
                    filt = parts[6]
                elif workflow == 'CLC':
                    filt = sample_type

                # Improved FUNCOTATION parsing
                gene = ''
                variant_class = ''
                if "FUNCOTATION=" in info:
                    try:
                        funcotation_str = info.split("FUNCOTATION=")[1]
                        funcotation_str = funcotation_str.strip('[]')
                        funcotation_entries = funcotation_str.split('],[')
                        for entry in funcotation_entries:
                            entry = entry.strip('[]')
                            funcotation_values = entry.split('|')
                            if len(funcotation_values) > 5:
                                gene_candidate = funcotation_values[0].strip().upper()
                                variant_class_candidate = funcotation_values[5].strip().upper()
                                # Only assign if not empty
                                if gene_candidate and variant_class_candidate:
                                    gene = gene_candidate
                                    variant_class = variant_class_candidate
                                    break  # Use the first valid annotation
                    except Exception as e:
                        print(f"Warning: Failed to parse FUNCOTATION in line with error: {e}")

                for sample in samples:
                    data.append([
                        chrom,
                        int(pos),
                        ref,
                        alt.split(',')[0],
                        sample,
                        workflow,
                        filt,
                        gene,
                        variant_class
                    ])


    return pd.DataFrame(data, columns=[
        'chr', 'pos', 'ref', 'alt',
        'sample_nr', 'workflow', 'filter_status',
        'gene', 'variant_class'
    ])



# Paths to the VCF files
gatk_path = "/mnt/scratch/nazilak/GATK/Variant_calling/homemade_PON/Liftover_PON/"
clc_annotated_path = "/mnt/scratch/nazilak/CLC_vcf/Funcotator/"
clc_low_freq_path = "/mnt/scratch/nazilak/CLC_vcf/Somatic_variants_with_low_frequency_in_germline/"

# Function to process all files
def process_files(gatk_path, clc_annotated_path, clc_low_freq_path):
    all_data = []

    for vcf_file in glob.glob(os.path.join(gatk_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "GATK")
        df['sample_nr'] = sample_name
        all_data.append(df)

    for vcf_file in glob.glob(os.path.join(clc_annotated_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "annotated")
        df['sample_nr'] = sample_name
        all_data.append(df)

    for vcf_file in glob.glob(os.path.join(clc_low_freq_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "low_frequency")
        df['sample_nr'] = sample_name
        all_data.append(df)


    combined = pd.concat(all_data, ignore_index=True)

    duplicates_within_sample = combined[combined.duplicated(subset=[
        'chr', 'pos', 'ref', 'alt', 'sample_nr', 'workflow', 'filter_status'
    ], keep=False)]
    print("Duplicate rows found within the same sample:")
    print(duplicates_within_sample)

    return combined

# Process and save
combined_data = process_files(gatk_path, clc_annotated_path, clc_low_freq_path)
combined_data.to_csv("/mnt/scratch/nazilak/Results/variant_comparisonHP.csv", index=False)
print("Combined data saved to /mnt/scratch/nazilak/Results/variant_comparisonHP.csv")


################################## SUPPLEMENTARY 3 ####################################
# ---------------------------------------------
# Summarize and plot GATK filter information
gatk_filter_summary_path = "/mnt/scratch/nazilak/Results/GATK_filters_overviewHP.txt"
plot_path = "/mnt/scratch/nazilak/Results/GATK_filters_overviewHP.png"

gatk_data = combined_data[combined_data['workflow'] == 'GATK']
filter_counts_with_samples = {}

for _, row in gatk_data.iterrows():
    filters = row['filter_status'].split(';')
    sample_nr = row['sample_nr']
    for filter_item in filters:
        if filter_item not in filter_counts_with_samples:
            filter_counts_with_samples[filter_item] = {'count': 0, 'samples': set()}
        filter_counts_with_samples[filter_item]['count'] += 1
        filter_counts_with_samples[filter_item]['samples'].add(sample_nr)

with open(gatk_filter_summary_path, "w") as filter_file:
    filter_file.write("GATK Individual Filtering Summary with Samples:\n")
    filter_file.write("===============================================\n")
    for filter_item, info in sorted(filter_counts_with_samples.items(), key=lambda x: x[1]['count'], reverse=True):
        samples_list = ', '.join(sorted(info['samples']))
        filter_file.write(f"{filter_item}: {info['count']} occurrences\n")
        filter_file.write(f"Found in samples: {samples_list}\n\n")

print(f"GATK individual filtering summary with samples saved to {gatk_filter_summary_path}")

# Prepare data
filters = []
counts = []
for filter_item, info in sorted(filter_counts_with_samples.items(), key=lambda x: x[1]['count'], reverse=True):
    filters.append(filter_item)
    counts.append(info['count'])

# Create a more compact and readable bar chart
# Create a taller and narrower bar chart
plt.figure(figsize=(7, 6))  # Reduced width, increased height

bars = plt.bar(filters, counts, color='steelblue')

plt.xlabel('Filtering Criteria', fontsize=16, fontweight='bold')
plt.ylabel('Number of Occurrences', fontsize=16, fontweight='bold')

plt.xticks(rotation=80, fontsize=13, ha='right')
plt.yticks(fontsize=13)

plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout(pad=1.5)
plt.savefig(plot_path, dpi=150)

print(f"Bar chart saved to {plot_path}")

################################ SUPPLEMENTARY 4.A ##################################
# ---------------------------------------------
# Create heatmap of GATK filters per sample (relative)

heatmap_path = "/mnt/scratch/nazilak/Results/gatk_filter_heatmap_relative_ALL_HP.png"

df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparisonHP.csv")
df_gatk = df[df['workflow'] == 'GATK']
df_expanded = df_gatk.assign(filter_status=df_gatk['filter_status'].str.split(';')).explode('filter_status')
heatmap_data = df_expanded.groupby(['filter_status', 'sample_nr']).size().unstack(fill_value=0)
total_variants_per_sample = df_gatk.groupby('sample_nr').size()
heatmap_data_relative = (heatmap_data / total_variants_per_sample) * 100
heatmap_data_relative = heatmap_data_relative.sort_index()

# Make the figure less tall and wider
plt.figure(figsize=(10, 5))  # Adjust width and height as needed

sns.heatmap(
    heatmap_data_relative,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    annot=False,
    cbar_kws={'label': 'Percentage (%)'}
)

plt.xlabel("Sample Number", fontsize=14)
plt.ylabel("Filter Type", fontsize=14)

# Show only every 5th sample number on X axis for readability
step = 5
xticks = range(0, len(heatmap_data_relative.columns), step)
xticklabels = heatmap_data_relative.columns[::step]
plt.xticks(
    ticks=xticks,
    labels=xticklabels,
    rotation=90,
    fontsize=12
)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(heatmap_path, dpi=200)
plt.show()

print(f"Heatmap saved as '{heatmap_path}'")

################################ SUPPLEMENTARY 4.B ##################################
###########################################################
# Create heatmap of GATK filters per sample (relative) from only FFPE Blocks
# Define the list of sample numbers to include (total 47)
selected_samples = {
    1, 4, 5, 8, 9, 11, 13, 14, 17, 18, 19, 20,
    22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33,
    34, 36, 37, 38, 39, 42, 43, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 57, 58, 59, 60, 61, 62, 63
}
selected_samples = {str(s) for s in selected_samples}  # ensure matching string type

filtered_heatmap_path = "/mnt/scratch/nazilak/Results/gatk_filter_heatmap_relative_HP_FFPE.png"

# Load the data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparisonHP.csv")

# Filter for GATK workflow and selected sample numbers
df_gatk = df[(df['workflow'] == 'GATK') & (df['sample_nr'].astype(str).isin(selected_samples))]

# Expand filter_status into multiple rows
df_expanded = df_gatk.assign(filter_status=df_gatk['filter_status'].str.split(';')).explode('filter_status')

# Group and compute relative frequencies
heatmap_data = df_expanded.groupby(['filter_status', 'sample_nr']).size().unstack(fill_value=0)
total_variants_per_sample = df_gatk.groupby('sample_nr').size()
heatmap_data_relative = (heatmap_data / total_variants_per_sample) * 100
heatmap_data_relative = heatmap_data_relative.sort_index()

# Make the figure less tall and wider
plt.figure(figsize=(9, 5))  # Adjust width and height as needed

sns.heatmap(
    heatmap_data_relative,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    annot=False,
    cbar_kws={'label': 'Percentage (%)'}
)

plt.xlabel("Sample Number", fontsize=14)
plt.ylabel("Filter Type", fontsize=14)

# Show only every 5th sample number on X axis for readability
step = 5
xticks = range(0, len(heatmap_data_relative.columns), step)
xticklabels = heatmap_data_relative.columns[::step]
plt.xticks(
    ticks=xticks,
    labels=xticklabels,
    rotation=90,
    fontsize=12
)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(filtered_heatmap_path, dpi=200)
plt.show()

print(f"Heatmap saved as '{filtered_heatmap_path}'")

################################ SUPPLEMENTARY 4.C ##################################
# Create heatmap of GATK filters per sample (relative) from only SMEARS
#Total 16 samples
selected_samples = {
    2, 3, 6, 7, 10, 12, 15, 16, 21, 28, 
    35, 40, 41, 44, 55, 56
}
selected_samples = {str(s) for s in selected_samples}  # ensure matching string type

filtered_heatmap_path = "/mnt/scratch/nazilak/Results/gatk_filter_heatmap_relative_HP_SMEAR.png"

# Load the data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparisonHP.csv")

# Filter for GATK workflow and selected sample numbers
df_gatk = df[(df['workflow'] == 'GATK') & (df['sample_nr'].astype(str).isin(selected_samples))]

# Expand filter_status into multiple rows
df_expanded = df_gatk.assign(filter_status=df_gatk['filter_status'].str.split(';')).explode('filter_status')

# Group and compute relative frequencies
heatmap_data = df_expanded.groupby(['filter_status', 'sample_nr']).size().unstack(fill_value=0)
total_variants_per_sample = df_gatk.groupby('sample_nr').size()
heatmap_data_relative = (heatmap_data / total_variants_per_sample) * 100
heatmap_data_relative = heatmap_data_relative.sort_index()

# Make the figure less tall and wider
plt.figure(figsize=(8, 5))  # Adjust width and height as needed

sns.heatmap(
    heatmap_data_relative,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    annot=False,
    cbar_kws={'label': 'Percentage (%)'}
)

plt.xlabel("Sample Number", fontsize=14)
plt.ylabel("Filter Type", fontsize=14)

# Show only every 5th sample number on X axis for readability
step = 5
xticks = range(0, len(heatmap_data_relative.columns), step)
xticklabels = heatmap_data_relative.columns[::step]
plt.xticks(
    ticks=xticks,
    labels=xticklabels,
    rotation=90,
    fontsize=12
)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(filtered_heatmap_path, dpi=200)
plt.show()
print(f"Filtered heatmap saved as '{filtered_heatmap_path}'")