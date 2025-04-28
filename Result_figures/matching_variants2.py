import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from contextlib import redirect_stdout

# Define the output file for terminal logs
log_file_path = "/mnt/scratch/nazilak/Results/Matching_variants/matching_variants2_log.txt"

# Load the CSV file containing the variant data
df = pd.read_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv")

# Filter the data to focus on GATK (PASS) and CLC (annotated) variants
gatk = df[df['workflow'] == 'GATK'].copy()
clc = df[df['workflow'] == 'CLC'].copy()

# Now it's safe to modify
gatk['chr'] = gatk['chr'].str.replace('chr', '', regex=False)
clc['chr'] = clc['chr'].str.replace('chr', '', regex=False)
gatk['pos'] = gatk['pos'].astype(int)
clc['pos'] = clc['pos'].astype(int)
# Remove duplicates (if any)
gatk = gatk.drop_duplicates(subset=['chr', 'pos', 'ref', 'alt', 'sample_nr'])
clc = clc.drop_duplicates(subset=['chr', 'pos', 'ref', 'alt', 'sample_nr'])

# Initialize variables to store results for the stacked bar plot
gatk_unique_counts = []  # Number of unique variants from GATK per sample
clc_unique_counts = []   # Number of unique variants from CLC per sample
match_counts = []        # Number of matching variants per sample

samples = range(1, 64)  # Sample numbers from 1 to 63

# Redirect stdout to both terminal and log file
with open(log_file_path, "w") as log_file:
    with redirect_stdout(log_file):
        # Check each sample for matching and unique variants
        all_matches = pd.DataFrame()  # To store all matches across samples

        for sample in samples:
            gatk_sample = gatk[gatk['sample_nr'] == sample]
            clc_sample = clc[clc['sample_nr'] == sample]
            
            # Merge based on matching chr, pos, ref, and alt
            matches = pd.merge(
                gatk_sample,
                clc_sample,
                on=['chr', 'pos', 'ref', 'alt'],
                how='inner'
            )
            
            # Count matches
            match_count = len(matches)
            match_counts.append(match_count)
            
            # Calculate unique GATK variants
            unique_gatk = len(gatk_sample) - match_count
            
            # Calculate unique CLC variants
            unique_clc = len(clc_sample) - match_count
            
            # Append the counts for unique variants
            gatk_unique_counts.append(unique_gatk)
            clc_unique_counts.append(unique_clc)

            # Append matches to the global DataFrame
            all_matches = pd.concat([all_matches, matches])

            # Print matches for this sample (if any)
            if not matches.empty:
                print(f"Sample {sample} matches:")
                print(matches)

        # Stacked bar plot: Three parts per sample
        plt.figure(figsize=(12, 10))

        # Stacked bar plot: Three parts per sample
        plt.bar(samples, gatk_unique_counts, color='lightblue', label='Unique GATK Variants')
        plt.bar(samples, clc_unique_counts, bottom=gatk_unique_counts, color='lightgreen', label='Unique CLC Variants')
        plt.bar(samples, match_counts, bottom=[i+j for i,j in zip(gatk_unique_counts, clc_unique_counts)], color='salmon', label='Matching Variants')

        # Labels and title
        plt.xlabel("Sample Number")
        plt.ylabel("Number of Variants")
        plt.title("Stacked Bar Plot of Matching and Unique Variants Per Sample (GATK (ALL) CLC (ALL))")

        # Add legend
        plt.legend()

        # Rotate x-axis labels for clarity
        plt.xticks(samples, rotation=90)
        plt.tight_layout()

        # Save the stacked bar plot
        plot_path = "/mnt/scratch/nazilak/Results/Matching_variants/stacked_matching_variants2_per_sample.png"
        plt.savefig(plot_path)

        # Print plot save location
        print(f"Stacked bar plot saved to {plot_path}")

        # Save match counts to CSV (optional)
        csv_path = "/mnt/scratch/nazilak/Results/Matching_variants/matching_variants2_per_sample.csv"
        output_df = pd.DataFrame({'sample_nr': samples, 'matching_variants': match_counts})
        output_df.to_csv(csv_path, index=False)

        # Print CSV save location
        print(f"Match counts saved to {csv_path}")

        # Venn Diagram: Total Variants vs Matching Variants
        total_gatk_variants = len(gatk)
        total_clc_variants = len(clc)
        total_matching_variants = len(all_matches)

        unique_gatk_variants = total_gatk_variants - total_matching_variants
        unique_clc_variants = total_clc_variants - total_matching_variants

        plt.figure(figsize=(8, 8))
        venn = venn2(
            subsets=(unique_gatk_variants, unique_clc_variants, total_matching_variants),
            set_labels=('GATK Total Variants', 'CLC Total Variants'),
            set_colors=('lightblue', 'lightgreen')  # Initial fill for the sets
        )

        # Customize intersection (matching variants) color to salmon
        venn.get_patch_by_id('11').set_color('salmon')  # '11' is the overlapping area
        venn.get_patch_by_id('10').set_color('lightblue')  # Unique GATK
        venn.get_patch_by_id('01').set_color('lightgreen')  # Unique CLC
        plt.title("Venn Diagram: GATK (ALL) vs CLC Variants (ALL)")

        venn_plot_path = "/mnt/scratch/nazilak/Results/Matching_variants/venn_diagram2_gatk_clc.png"
        plt.savefig(venn_plot_path)

        # Print Venn diagram save location
        print(f"Venn diagram saved to {venn_plot_path}")

        # Load the CSV file containing the matching variants per sample
        df_matching = pd.read_csv("/mnt/scratch/nazilak/Results/Matching_variants/matching_variants2_per_sample.csv")

        # Create a bar plot of matching variants per sample
        plt.figure(figsize=(12, 5))
        plt.bar(df_matching['sample_nr'], df_matching['matching_variants'], color='lightblue')
        plt.xlabel("Sample Number")
        plt.ylabel("Number of Common Variants")
        plt.title("Common Variants Between GATK (ALL) and CLC (ALL) Per Sample")
        plt.xticks(samples, rotation=90)
        plt.tight_layout()

        # Save the bar plot for common variants
        common_variants_plot_path = "/mnt/scratch/nazilak/Results/Matching_variants/common_variants_per_sample2.png"
        plt.savefig(common_variants_plot_path)

        # Print bar plot save location
        print(f"Bar plot for common variants saved to {common_variants_plot_path}")
