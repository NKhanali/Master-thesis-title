import os
import re
import csv
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

# Updated Function to process each VCF file
def process_vcf_file(vcf_file_path, sample_number, csv_writer):
    with open(vcf_file_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')

            info_field = columns[7]

            # Extract FUNCOTATION from INFO column
            match = re.search(r'FUNCOTATION=\[([^\]]+)\]', info_field)
            if match:
                funcotation = match.group(1)
                func_parts = funcotation.split('|')

                # Ensure the expected number of fields exist
                if len(func_parts) >= 10:
                    selected_fields = [
                        func_parts[0],  # Gene name
                        func_parts[5],  # Variant Classification
                        func_parts[7],  # Variant Type
                    ]

                    row = [sample_number] + selected_fields
                    csv_writer.writerow(row)

# New input/output paths
vcf_directory = '/mnt/scratch/nazilak/CLC_vcf/Funcotator'  # Updated path
output_csv_file = '/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv'

# Write the CSV header and process VCFs
with open(output_csv_file, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)

    header = ['Sample_Nr', 'Gene name', 'Variant Classification', 'Variant Type']
    csv_writer.writerow(header)

    for filename in os.listdir(vcf_directory):
        if filename.endswith('_annotatedCLC.vcf'):
            sample_number = filename.split('_')[0]  # extract number 1-63
            vcf_file_path = os.path.join(vcf_directory, filename)
            process_vcf_file(vcf_file_path, sample_number, csv_writer)

print(f'Done! Output saved to {output_csv_file}')

################################################################
import pandas as pd
import matplotlib.pyplot as plt

# Load the filtered CSV file
csv_file = '/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv'
df = pd.read_csv(csv_file)

# Count occurrences of each Variant Classification
classification_counts = df['Variant Classification'].value_counts()

# Plot the bar chart
plt.figure(figsize=(12, 6))
classification_counts.plot(kind='bar', color='skyblue', edgecolor='black')

# Add labels and title
plt.xlabel('Variant Classification', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title('Count of Variant Classifications Across Samples CLC', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

os.makedirs('/mnt/scratch/nazilak/Results/Protein_change/', exist_ok=True)
plt.savefig('/mnt/scratch/nazilak/Results/Protein_change/barplot_variant_classificationCLC.png')
# Show plot
plt.show()
print(f'Done! barplot_variant_classification.png')

################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm


# Read the CSV file
def create_variant_heatmap(csv_path):
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Divide similar variant classification into the same category
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
        "COULD_NOT_DETERMINE": "UNKNOWN", 
        "THREE_PRIME_UTR": "NON_CODING",
        "DE_NOVO_START_OUT_FRAME": "FRAME_SHIFT",
        "START_CODON_INS": "FRAME_SHIFT",
        "START_CODON_SNP": "FRAME_SHIFT",
    }
    # Apply the classification map
    df["Merged Classification"] = df["Variant Classification"].map(classification_map)

    # Define a discrete rank for the merged categories
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
    
    # For each sample and gene, select the variant with the highest impact (lowest rank)
    selected_variants = df.loc[df.groupby(["Sample_Nr", "Gene name"])["Impact Rank"].idxmin()]
    
    # Create a pivot table using the impact ranks
    pivot_data = selected_variants.pivot(index="Gene name", columns="Sample_Nr", values="Impact Rank")
    
    # Get unique gene names from the original DataFrame
    unique_genes = df['Gene name'].unique()
    
    # Reindex the pivot data to include all unique genes
    pivot_data = pivot_data.reindex(unique_genes)
    
    # Sort genes by the number of samples with a variant (descending), most frequent at the bottom
    variant_counts = pivot_data.notna().sum(axis=1)
    pivot_data = pivot_data.loc[variant_counts.sort_values(ascending=True).index]

    # Ensure all sample numbers are included in the heatmap
    sample_numbers = sorted(df["Sample_Nr"].unique())
    
    # Reindex the columns to include all sample numbers
    pivot_data = pivot_data.reindex(columns=sample_numbers)

    # Convert column names to numeric sample numbers for proper sorting
    sample_numbers = df["Sample_Nr"].astype(str).str.extract(r"^(\d+)")[0].dropna().astype(int).unique()
    sample_numbers = sorted(sample_numbers)

    # Rename pivot columns to match numeric sample numbers
    pivot_data.columns = pivot_data.columns.astype(str).str.extract(r"^(\d+)")[0].astype(float).astype('Int64')
    pivot_data = pivot_data.reindex(columns=sample_numbers)

        
    # Create the heatmap
    plt.figure(figsize=(20, len(pivot_data.index) * 0.2))  # adjust 0.3 as needed

    # Generate a color palette with one color per unique rank
    num_categories = len(impact_ranking)
    palette = sns.color_palette("RdYlBu", n_colors=num_categories)

    # Create a discrete colormap
    cmap = ListedColormap(palette)

    # Define boundaries between categories
    bounds = list(range(1, num_categories + 2))  # +2 because bounds are edges, e.g. [1,2,3,...,11]
    norm = BoundaryNorm(bounds, cmap.N)

    # Create the heatmap using discrete colormap and norm
    ax = sns.heatmap(pivot_data, cmap=cmap, norm=norm, cbar_kws={'label': 'Variant Classification'})

    # Get the colorbar
    cbar = ax.collections[0].colorbar
    
    # Create a list of ranks and their corresponding variant classifications
    rank_to_variant = {rank: variant for variant, rank in impact_ranking.items()}
    
    # Show all defined variant types in the colorbar
    all_ranks = sorted(impact_ranking.values())

    # Set new tick positions for all ranks
    cbar.set_ticks(all_ranks)

    # Set tick labels to variant classifications, including UNKNOWN
    cbar.set_ticklabels([rank_to_variant[rank] for rank in all_ranks])

    # Label the colorbar
    cbar.set_label("Variant Classification")
    
    # Final formatting
    plt.title("Heatmap of Most Impactful Variant Classifications GATK")
    plt.xlabel("Sample Number")
    plt.ylabel("Gene Name")
    plt.xticks(fontsize=8)
    plt.tight_layout()

    output_dir = "/mnt/scratch/nazilak/Results/Protein_change/"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "heatmapCLC.svg"), format='svg', dpi=300, bbox_inches='tight')
    plt.show()

    print('Done! heatmap.svg')

# Example usage with your CSV file
create_variant_heatmap("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv")



######################################################
import csv

# Use a set to store unique gene names
unique_genes = set()

# Open your CSV file
with open("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv", mode="r", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        unique_genes.add(row["Gene name"])

# Write the unique gene names to list.txt
with open("/mnt/scratch/nazilak/Results/Protein_change/gene_listCLC.txt", mode="w") as outfile:
    for gene in sorted(unique_genes):
        outfile.write(gene + "\n")
print("Done! gene_list.txt")






