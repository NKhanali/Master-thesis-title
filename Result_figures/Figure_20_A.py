import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from collections import Counter


file1 = "/mnt/scratch/nazilak/Results/variant_comparison.csv"
file2 = "/mnt/scratch/nazilak/Results/variant_comparisonHP.csv"


df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

df1_gatk = df1[df1["workflow"] == "GATK"]
df2_gatk = df2[df2["workflow"] == "GATK"]


variants1 = df1_gatk.apply(lambda row: f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}", axis=1).tolist()
variants2 = df2_gatk.apply(lambda row: f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}", axis=1).tolist()


counter1 = Counter(variants1)
counter2 = Counter(variants2)


multiset1 = list(counter1.elements())
multiset2 = list(counter2.elements())


both = Counter(multiset1) & Counter(multiset2) 
only1 = Counter(multiset1) - Counter(multiset2)
only2 = Counter(multiset2) - Counter(multiset1)


both_count = sum(both.values())
only1_count = sum(only1.values())
only2_count = sum(only2.values())


print(f"Total in file 1 (GATK): {len(multiset1)}")
print(f"Total in file 2 (GATK): {len(multiset2)}")
print(f"Shared calls         : {both_count}")
print(f"Unique to file 1     : {only1_count}")
print(f"Unique to file 2     : {only2_count}")


plt.figure(figsize=(8, 8))
v = venn2(subsets=(only1_count, only2_count, both_count), set_labels=('', ''))


v.get_label_by_id('10').set_text(str(only1_count))
v.get_label_by_id('01').set_text(str(only2_count))
v.get_label_by_id('11').set_text(str(both_count))


plt.text(-0.7, 0.4, 'GATK (Default)', fontsize=14, ha='center')
plt.text(0.7, 0.4, 'GATK (HP)', fontsize=14, ha='center')

plt.title("All GATK Variant Calls (Including Duplicates)")
plt.tight_layout()
plt.savefig("/mnt/scratch/nazilak/Results/PON_vs_HP_ALL.png")
plt.show()



