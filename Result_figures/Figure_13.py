import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Paths to your gene list text files
#these two files were outputs from the supplementary 11 and 12 script
file1 = '/mnt/scratch/nazilak/Results/Protein_change/top20_CLC.txt'
file2 = '/mnt/scratch/nazilak/Results/Protein_change/top20_GATK.txt'

# Function to read gene names from a text file (one gene per line)
def read_genes(file_path):
    with open(file_path, 'r') as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes

# Read gene sets
genes1 = read_genes(file1)
genes2 = read_genes(file2)

# Create Venn diagram
venn = venn2([genes1, genes2], set_labels=('CLC', 'GATK'))

# Set colors
venn.get_patch_by_id('10').set_color('lightgreen')  # Only in Dataset 1
venn.get_patch_by_id('01').set_color('lightblue')   # Only in Dataset 2
venn.get_patch_by_id('11').set_color('salmon')      # Overlap

# Function to format gene names for labels (multi-line)
def format_label(gene_set):
    return '\n'.join(sorted(gene_set)) if gene_set else ''

# Set labels inside the diagram
venn.get_label_by_id('10').set_text(format_label(genes1 - genes2))
venn.get_label_by_id('01').set_text(format_label(genes2 - genes1))
venn.get_label_by_id('11').set_text(format_label(genes1 & genes2))

# Adjust layout and save figure
plt.tight_layout()
plt.savefig('/mnt/scratch/nazilak/Results/Protein_change/top20_match.png', dpi=150)
plt.show()
