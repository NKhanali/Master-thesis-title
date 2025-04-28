
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns


#GATK: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  1_N     1_T
#CLC:  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Annotated_somatic_variants-B025,_B099

# Define a function that analyzes the different VCF files and their path.
# Has 3 parameters: file_path (path to the VCF files), workflow  (GATK or CLC), and sample_type (used for CLC to separate annotated and lowfreq)
def vcf(file_path, workflow, sample_type=None):
    data = []  # A list for each sample for every row
    samples = []  # Holds sample IDs found in the VCF

    # Opens the file and "r"eads it
    with open(file_path, 'r') as f:
        for line in f:  # Read all the lines in the VCF file
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = list(set([sample.split('_')[0] for sample in parts[9:]])) #extract unique numeric sample IDs (ex 1 from 1_N and 1_T), remove duplicates

            elif not line.startswith('#'):  # For lines that do not start with #, these are variant data rows.
                parts = line.strip().split('\t')
                chrom = parts[0]
                if workflow == 'CLC' and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom #added chrom to the CLC contigs
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]

                if workflow == 'GATK':
                    filt = parts[6]  # The filters from GATK 
                elif workflow == 'CLC':
                    filt = sample_type  #whether the file name is annotated or low freq 

                for sample in samples:  # Loop creates a row per variant per sample. Each row contains:
                    data.append([
                        chrom,
                        int(pos),
                        ref,
                        alt.split(',')[0],
                        sample,
                        workflow,
                        filt
                    ])

    return pd.DataFrame(data, columns=[
        'chr', 'pos', 'ref', 'alt',
        'sample_nr', 'workflow', 'filter_status'
    ])

# Paths to the VCF files
gatk_path = "/mnt/scratch/nazilak/GATK/Variant_calling/Liftover/"
clc_annotated_path = "/mnt/scratch/nazilak/CLC_vcf/Annotated_somatic_variants/"
clc_low_freq_path = "/mnt/scratch/nazilak/CLC_vcf/Somatic_variants_with_low_frequency_in_germline/"


# This function processes all VCF files in three directories (GATK, CLC Annotated, CLC Low Frequency).
def process_files(gatk_path, clc_annotated_path, clc_low_freq_path):
    all_data = []

    # Process GATK files
    for vcf_file in glob.glob(os.path.join(gatk_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]  # Extracts the sample name (e.g., from 1_liftover.vcf, gets 1).
        df = vcf(vcf_file, "GATK")  # Uses the vcf function from above. Calls the vcf() function with workflow set to 'GATK'
        df['sample_nr'] = sample_name  # Adds the sample name (1, 2, 3 etc)
        all_data.append(df)  # Appends the dataframe to the full list

    # Process CLC annotated files
    for vcf_file in glob.glob(os.path.join(clc_annotated_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "annotated")
        df['sample_nr'] = sample_name
        all_data.append(df)

    # Process CLC low frequency files
    for vcf_file in glob.glob(os.path.join(clc_low_freq_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "low_frequency")
        df['sample_nr'] = sample_name
        all_data.append(df)

    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)

    #chck for duplicates
    duplicates_within_sample = combined[combined.duplicated(subset=['chr', 'pos', 'ref', 'alt', 'sample_nr', 'workflow', 'filter_status'], keep=False)]
    print("Duplicate rows found within the same sample:")
    print(duplicates_within_sample)

    return combined


# Process all files and combine the data
combined_data = process_files(gatk_path, clc_annotated_path, clc_low_freq_path)

# Save the combined data to a CSV file
combined_data.to_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv", index=False)

print("Combined data saved to /mnt/scratch/nazilak/Results/variant_comparison.csv")


################################################################################################

gatk_filter_summary_path = "/mnt/scratch/nazilak/Results/GATK_filters_overview.txt"
plot_path = "/mnt/scratch/nazilak/Results/GATK_filters_overview.png"


gatk_data = combined_data[combined_data['workflow'] == 'GATK'] #only gatk variants

# Track filters and associated samples
filter_counts_with_samples = {}
for _, row in gatk_data.iterrows():
    filters = row['filter_status'].split(';')  # Split by ';' so each filter is individual
    sample_nr = row['sample_nr'] #include the sample nr of each filter to calculate which samples each filter were on
    
    for filter_item in filters:
        if filter_item not in filter_counts_with_samples: #doesnt 
            filter_counts_with_samples[filter_item] = {'count': 0, 'samples': set()} 
            #'count': 0: A counter set to zero (to track how many times this filter has been applied).
            #'samples': set(): An empty set (to store unique sample numbers where this filter was applied).
        
        filter_counts_with_samples[filter_item]['count'] += 1
        filter_counts_with_samples[filter_item]['samples'].add(sample_nr)

# **Step 1: Write Individual Filter Summary to Text File**
with open(gatk_filter_summary_path, "w") as filter_file:
    filter_file.write("GATK Individual Filtering Summary with Samples:\n")
    filter_file.write("===============================================\n")
    for filter_item, info in sorted(filter_counts_with_samples.items(), key=lambda x: x[1]['count'], reverse=True):
        samples_list = ', '.join(sorted(info['samples']))  # Sort sample numbers for readability
        filter_file.write(f"{filter_item}: {info['count']} occurrences\n")
        filter_file.write(f"Found in samples: {samples_list}\n\n")

print(f"GATK individual filtering summary with samples saved to {gatk_filter_summary_path}")

# **Step 2: Prepare Data for Bar Chart**
filters = []
counts = []
for filter_item, info in sorted(filter_counts_with_samples.items(), key=lambda x: x[1]['count'], reverse=True):
    filters.append(filter_item)
    counts.append(info['count'])

# Create a bar chart
plt.figure(figsize=(12, 6))
plt.bar(filters, counts, color='skyblue')
plt.xlabel('Filtering Criteria')
plt.ylabel('Number of Occurrences')
plt.title('GATK Filtering Criteria and Their Counts')
plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
plt.tight_layout()

# Save the figure
plt.savefig(plot_path)
print(f"Bar chart saved to {plot_path}")

################################################################################################
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the combined data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")

# Step 1: Filter for GATK workflow
df_gatk = df[df['workflow'] == 'GATK']
# Step 2: Split the semicolon-separated filters into individual rows
df_expanded = df_gatk.assign(filter_status=df_gatk['filter_status'].str.split(';')).explode('filter_status')

# Step 3: Create a pivot table counting total occurrences of each filter in each sample
heatmap_data = df_expanded.groupby(['filter_status', 'sample_nr']).size().unstack(fill_value=0)

# Step 4: Calculate the total number of variants per sample
# This is the count of non-filtered variants per sample
total_variants_per_sample = df_gatk.groupby('sample_nr').size()

# Step 5: Normalize values to percentages relative to total variants per sample
# Normalize the heatmap data by dividing each value by the total number of variants in the corresponding sample
heatmap_data_relative = (heatmap_data / total_variants_per_sample) * 100

# Optional: Sort for nicer visual appearance
heatmap_data_relative = heatmap_data_relative.sort_index()

# Plot the heatmap with relative values
plt.figure(figsize=(27, 10))
sns.heatmap(
    heatmap_data_relative,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    annot=True, # Show percentages in each cell
    fmt='.1f', # Format annotations to one decimal place
    cbar_kws={'label': 'Percentage (%)'}
)
plt.title("GATK Filter Occurrences per Sample (Relative)", fontsize=16)
plt.xlabel("Sample Number")
plt.ylabel("Filter Type")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Save the figure to a file
heatmap_path = "/mnt/scratch/nazilak/Results/gatk_filter_heatmap_relative_ALL.png"
plt.savefig(heatmap_path, dpi=300) # Save as PNG

print(f"Heatmap saved as '{heatmap_path}'")
