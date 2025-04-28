import os
import re
import csv
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

# Function to process each VCF file
def process_vcf_file(vcf_file_path, sample_number, csv_writer):
    with open(vcf_file_path, 'r') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')

            if columns[6] != 'PASS': #only the passed variants
                continue 

            info_field = columns[7]

            # Extract FUNCOTATION from INFO column
            match = re.search(r'FUNCOTATION=\[([^\]]+)\]', info_field)
            if match:
                funcotation = match.group(1)
                func_parts = funcotation.split('|')

                # Make sure we have enough fields
                if len(func_parts) >= 10:
                    # Extract only selected FUNCOTATION fields
                    selected_fields = [
                        func_parts[0],  # FUNCOTATION_1
                        func_parts[5],  # FUNCOTATION_6
                        func_parts[7],  # FUNCOTATION_8
                    ]

                    # Write row with sample number and selected fields
                    row = [sample_number] + selected_fields
                    csv_writer.writerow(row)
    


# Input and output paths
vcf_directory = '/mnt/scratch/nazilak/GATK/Variant_calling/Liftover'  # Update as needed
output_csv_file = '/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_output.csv'


# Open CSV for writing
with open(output_csv_file, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Header row with selected columns
    header = [
        'Sample_Nr',
        'Gene name',       # FUNCOTATION_1
        'Variant Classification', # FUNCOTATION_6cd ..
        'Variant Type' #FUNCOTATION_8

    ]
    csv_writer.writerow(header)

    # Loop through VCFs
    for filename in os.listdir(vcf_directory):
        if filename.endswith('_liftover.vcf'):
            sample_number = filename.split('_')[0]  # Extract number from 1-63
            vcf_file_path = os.path.join(vcf_directory, filename)
            #print(f'Processing: {filename}')
            process_vcf_file(vcf_file_path, sample_number, csv_writer)
print(f'Done! Output saved to {output_csv_file}')

################################################################
import pandas as pd
import matplotlib.pyplot as plt

# Load the filtered CSV file
csv_file = '/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_output.csv'
df = pd.read_csv(csv_file)

# Count occurrences of each Variant Classification
classification_counts = df['Variant Classification'].value_counts()

# Plot the bar chart
plt.figure(figsize=(12, 6))
classification_counts.plot(kind='bar', color='skyblue', edgecolor='black')

# Add labels and title
plt.xlabel('Variant Classification', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title('Count of Variant Classifications Across Samples', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

os.makedirs('/mnt/scratch/nazilak/Results/Protein_change/', exist_ok=True)
plt.savefig('/mnt/scratch/nazilak/Results/Protein_change/barplot_variant_classification.png')
# Show plot
plt.show()
print(f'Done! barplot_variant_classification.png')

################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm


def create_variant_heatmap(csv_path):
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Map similar variant classifications into broader categories
    classification_map = {
        "MISSENSE": "MISSENSE",
        "FRAME_SHIFT_DEL": "FRAME_SHIFT",
        "FRAME_SHIFT_INS": "FRAME_SHIFT",
        "SILENT": "SILENT",
        "SPLICE_SITE": "SPLICE_SITE",
        "INTRON": "NON_CODING",
        "NONSENSE": "NONSENSE",
        "IN_FRAME_DEL": "IN_FRAME",
        "IN_FRAME_INS": "IN_FRAME",
        "FIVE_PRIME_UTR": "NON_CODING",
        "COULD_NOT_DETERMINE": "UNKNOWN"
    }
    df["Merged Classification"] = df["Variant Classification"].map(classification_map)

    # Ranking impact of each classification
    impact_ranking = {
        "NONSENSE": 1,
        "FRAME_SHIFT": 2,
        "SPLICE_SITE": 3,
        "MISSENSE": 4,
        "IN_FRAME": 5,
        "SILENT": 6,
        "NON_CODING": 7,
        "UNKNOWN": 8
    }
    df["Impact Rank"] = df["Merged Classification"].map(impact_ranking)

    # Select the most impactful variant per gene per sample
    selected_variants = df.loc[df.groupby(["Sample_Nr", "Gene name"])["Impact Rank"].idxmin()]

    # Pivot table
    pivot_data = selected_variants.pivot(index="Gene name", columns="Sample_Nr", values="Impact Rank")

    # Reindex gene order based on how frequently they're affected
    unique_genes = df["Gene name"].unique()
    pivot_data = pivot_data.reindex(unique_genes)
    variant_counts = pivot_data.notna().sum(axis=1)
    pivot_data = pivot_data.loc[variant_counts.sort_values(ascending=True).index]

    # --- Fix the x-axis for compactness ---
    # Extract numeric sample numbers
    sample_numbers = df["Sample_Nr"].astype(str).str.extract(r"(\d+)")[0].astype(int).unique()
    sample_numbers = sorted(sample_numbers)

    # Rename and sort columns numerically
    pivot_data.columns = pivot_data.columns.astype(str).str.extract(r"(\d+)")[0].astype(float).astype('Int64')
    pivot_data = pivot_data.reindex(columns=sample_numbers)

    # Ensure all sample numbers (e.g. 1–63) are in the x-axis (this is because sample 55 and 51)
    full_sample_range = list(range(1, 64))  # 1 to 63 inclusive
    pivot_data = pivot_data.reindex(columns=full_sample_range)

    # Create the heatmap
    plt.figure(figsize=(20, len(pivot_data.index) * 0.2))

    num_categories = len(impact_ranking)
    palette = sns.color_palette("RdYlBu", n_colors=num_categories)
    cmap = ListedColormap(palette)
    bounds = list(range(1, num_categories + 2))
    norm = BoundaryNorm(bounds, cmap.N)

    ax = sns.heatmap(pivot_data, cmap=cmap, norm=norm, cbar_kws={'label': 'Variant Classification'})

    # Configure colorbar ticks and labels
    cbar = ax.collections[0].colorbar
    rank_to_variant = {rank: variant for variant, rank in impact_ranking.items()}
    all_ranks = sorted(impact_ranking.values())
    cbar.set_ticks(all_ranks)
    cbar.set_ticklabels([rank_to_variant[rank] for rank in all_ranks])
    cbar.set_label("Variant Classification")

    # Final formatting
    plt.title("Heatmap of Most Impactful Variant Classifications GATK")
    plt.xlabel("Sample Number")
    plt.ylabel("Gene Name")
    plt.xticks(fontsize=8)
    plt.tight_layout()

    output_dir = "/mnt/scratch/nazilak/Results/Protein_change/"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "heatmapGATK.svg"), format='svg', dpi=300, bbox_inches='tight')
    plt.show()

    print("Done! heatmapGATK.svg")

# Example usage
create_variant_heatmap("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_output.csv")



######################################################
import csv

# Use a set to store unique gene names
unique_genes = set()

# Open your CSV file
with open("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_output.csv", mode="r", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        unique_genes.add(row["Gene name"])

# Write the unique gene names to list.txt
with open("/mnt/scratch/nazilak/Results/Protein_change/gene_list.txt", mode="w") as outfile:
    for gene in sorted(unique_genes):
        outfile.write(gene + "\n")
print("Done! gene_list.txt")










