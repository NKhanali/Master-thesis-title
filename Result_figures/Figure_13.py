import matplotlib.pyplot as plt
from matplotlib_venn import venn2


file1 = '/mnt/scratch/nazilak/Results/Protein_change/top20_CLC.txt'
file2 = '/mnt/scratch/nazilak/Results/Protein_change/top20_GATK.txt'


def read_genes(file_path):
    with open(file_path, 'r') as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes


genes1 = read_genes(file1)
genes2 = read_genes(file2)

venn = venn2([genes1, genes2], set_labels=('CLC', 'GATK'))

venn.get_patch_by_id('10').set_color('lightgreen')  
venn.get_patch_by_id('01').set_color('lightblue')   
venn.get_patch_by_id('11').set_color('salmon')      


def format_label(gene_set):
    return '\n'.join(sorted(gene_set)) if gene_set else ''

venn.get_label_by_id('10').set_text(format_label(genes1 - genes2))
venn.get_label_by_id('01').set_text(format_label(genes2 - genes1))
venn.get_label_by_id('11').set_text(format_label(genes1 & genes2))

plt.tight_layout()
plt.savefig('/mnt/scratch/nazilak/Results/Protein_change/top20_match.png', dpi=150)
plt.show()
