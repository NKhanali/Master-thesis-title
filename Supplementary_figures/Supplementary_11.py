import os
import re
import csv
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Updated Function to process each VCF file
def process_vcf_file(vcf_file_path, sample_number, csv_writer):
    with open(vcf_file_path, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')

            info_field = columns[7]

            # Extract FUNCOTATION from INFO column
            match = re.search(r'FUNCOTATION=\[([^\]]+)\]', info_field)
            if match:
                funcotation = match.group(1)
                func_parts = funcotation.split('|')

                # Make sure we have enough fields
                if len(func_parts) >= 10:
                    selected_fields = [
                        func_parts[0],  # Gene name
                        func_parts[5],  # Variant Classification
                        func_parts[7],  # Variant Type
                    ]

                    # Write row with sample number and selected fields
                    row = [sample_number] + selected_fields
                    csv_writer.writerow(row)

# Input and output paths
vcf_directory = '/mnt/scratch/nazilak/CLC_vcf/Funcotator'  
output_csv_file = '/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv'

# Open CSV for writing
with open(output_csv_file, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Header row with selected columns
    header = ['Sample_Nr', 'Gene name', 'Variant Classification', 'Variant Type']
    csv_writer.writerow(header)

    # Loop through VCFs
    for filename in os.listdir(vcf_directory):
        if filename.endswith('_annotatedCLC.vcf'):
            sample_number = filename.split('_')[0]  # Extract number 1-63
            vcf_file_path = os.path.join(vcf_directory, filename)
            process_vcf_file(vcf_file_path, sample_number, csv_writer)
print(f'Done! Output saved to {output_csv_file}')

################################################################


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
    df["Merged Classification"] = df["Variant Classification"].map(classification_map)

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
    pivot_data = selected_variants.pivot(index="Gene name", columns="Sample_Nr", values="Impact Rank")
    
    unique_genes = df['Gene name'].unique()
    pivot_data = pivot_data.reindex(unique_genes)
    variant_counts = pivot_data.notna().sum(axis=1)
    pivot_data = pivot_data.loc[variant_counts.sort_values(ascending=True).index]

    # Ensure all sample numbers are included in the heatmap
    sample_numbers = df["Sample_Nr"].astype(str).str.extract(r"(\d+)")[0].astype(int).unique()
    sample_numbers = sorted(sample_numbers)
    pivot_data.columns = pivot_data.columns.astype(str).str.extract(r"(\d+)")[0].astype(float).astype('Int64')
    pivot_data = pivot_data.reindex(columns=sample_numbers)
    full_sample_range = list(range(1, 64))  # 1 to 63 inclusive
    pivot_data = pivot_data.reindex(columns=full_sample_range)

    # --- Begin Option 1: Custom colorbar ---
    fig, ax = plt.subplots(figsize=(20, len(pivot_data.index) * 0.2))

    num_categories = len(impact_ranking)
    palette = sns.color_palette("RdYlBu", n_colors=num_categories)
    cmap = ListedColormap(palette)
    bounds = list(range(1, num_categories + 2))
    norm = BoundaryNorm(bounds, cmap.N)

    # Draw heatmap without colorbar
    sns.heatmap(
        pivot_data,
        cmap=cmap,
        norm=norm,
        ax=ax,
        cbar=False
    )

    # Create custom colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(
        sm,
        ax=ax,
        orientation='vertical',
        fraction=0.03,  # smaller = thinner colorbar
        pad=0.02,
        aspect=120        # higher = thinner colorbar
    )

    # Configure colorbar ticks and labels
    rank_to_variant = {rank: variant for variant, rank in impact_ranking.items()}
    all_ranks = sorted(impact_ranking.values())
    cbar.set_ticks(all_ranks)
    cbar.set_ticklabels([rank_to_variant[rank] for rank in all_ranks])
    cbar.set_label("Variant Classification")

    # Final formatting
    ax.set_title("Heatmap of Most Impactful Variant Classifications CLC")
    ax.set_xlabel("Sample Number")
    ax.set_ylabel("Gene Name")
    ax.tick_params(axis='x', labelsize=8)
    plt.tight_layout()

    output_dir = "/mnt/scratch/nazilak/Results/Protein_change/"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "heatmapCLC.svg"), format='svg', dpi=300, bbox_inches='tight')
    plt.show()

    print('Done! heatmapCLC.svg')


# Example usage:
create_variant_heatmap("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv")



##########################################
# Load the filtered CSV file
df = pd.read_csv('/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_outputCLC.csv')

# Count variants per gene (each gene only once)
gene_counts = df['Gene name'].value_counts().head(20)

# Prepare output path in the same directory as the CSV
output_path = '/mnt/scratch/nazilak/Results/Protein_change/top20_CLC.txt'

# Write the top 20 genes to the text file (each gene only once)
with open(output_path, 'w') as f:
    for gene in gene_counts.index:
        f.write(f"{gene}\n")

print(f"Top 20 genes written to {output_path}")

