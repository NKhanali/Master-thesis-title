import os
import re
import csv
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

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

    # For each sample and gene, select the variant with the highest impact (lowest rank)
    selected_variants = df.loc[df.groupby(["Sample_Nr", "Gene name"])["Impact Rank"].idxmin()]
    pivot_data = selected_variants.pivot(index="Gene name", columns="Sample_Nr", values="Impact Rank")

    # Only show selected genes (the match between top 20 CLC and GATK genes)
    selected_genes = ["ARID1B", "CREBBP", "FOXL2", "KMT2C", "KMT2D", "LRP1B", "NOTCH1", "RB1", "ROS1", "TP53"]
    pivot_data = pivot_data.loc[selected_genes]

    variant_counts = pivot_data.notna().sum(axis=1)
    pivot_data = pivot_data.loc[variant_counts.sort_values(ascending=True).index]

    # --- Fix the x-axis for compactness ---
    # Extract numeric sample numbers
    sample_numbers = df["Sample_Nr"].astype(str).str.extract(r"(\d+)")[0].astype(int).unique()
    sample_numbers = sorted(sample_numbers)
    pivot_data.columns = pivot_data.columns.astype(str).str.extract(r"(\d+)")[0].astype(float).astype('Int64')
    pivot_data = pivot_data.reindex(columns=sample_numbers)
    full_sample_range = list(range(1, 64))  # 1 to 63 inclusive
    pivot_data = pivot_data.reindex(columns=full_sample_range)

    plt.figure(figsize=(7, 4))


   # --- Begin Option 1: Custom colorbar ---
    fig, ax = plt.subplots(figsize=(10,5))


    # Generate a color palette with one color per unique rank
    num_categories = len(impact_ranking)
    palette = sns.color_palette("RdYlBu", n_colors=num_categories)
    cmap = ListedColormap(palette)
    bounds = list(range(1, num_categories + 2))
    norm = BoundaryNorm(bounds, cmap.N)
    #ax = sns.heatmap(pivot_data, cmap=cmap, norm=norm, cbar_kws={'label': 'Variant Classification'})
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
        fraction=0.05,  # smaller = thinner colorbar
        pad=0.02,
        aspect=50        # higher = thinner colorbar
    )

    # Configure colorbar ticks and labels
    #cbar = ax.collections[0].colorbar
    rank_to_variant = {rank: variant for variant, rank in impact_ranking.items()}
    all_ranks = sorted(impact_ranking.values())
    cbar.set_ticks(all_ranks)
    cbar.set_ticklabels([rank_to_variant[rank] for rank in all_ranks])
    cbar.set_label("Variant Classification", fontsize=16)

    # Final formatting
    plt.xlabel("Sample Number", fontsize=16)
    plt.ylabel("Gene Name", fontsize=16)
    plt.xticks(fontsize=8)
    plt.tight_layout()

    output_dir = "/mnt/scratch/nazilak/Results/Protein_change/"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "heatmapGATK_top10.svg"), format='svg', dpi=150, bbox_inches='tight')
    plt.show()

    print("Done! heatmapGATK.svg")


create_variant_heatmap("/mnt/scratch/nazilak/Results/Protein_change/filtered_funcotation_output.csv")

