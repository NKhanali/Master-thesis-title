import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

ffpe_samples = {
    1, 4, 5, 8, 9, 11, 13, 14, 17, 18, 19, 20,
    22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33,
    34, 36, 37, 38, 39, 42, 43, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 57, 58, 59, 60, 61, 62, 63
}
ffpe_samples = {str(s) for s in ffpe_samples}

smear_samples = {
    2, 3, 6, 7, 10, 12, 15, 16, 21, 28, 
    35, 40, 41, 44, 55, 56
}
smear_samples = {str(s) for s in smear_samples}

csv_path = "/mnt/scratch/nazilak/Results/variant_comparison.csv"
df = pd.read_csv(csv_path)

df_gatk = df[df['workflow'] == 'GATK'].copy()


df_expanded = df_gatk.assign(filter_status=df_gatk['filter_status'].str.split(';')).explode('filter_status')

########################### FIGURE 10.A ##################################
def get_contamination_percentages(df_expanded, sample_set, group_name):
    data = []
    for sample in sample_set:
        sample_variants = df_expanded[df_expanded['sample_nr'].astype(str) == sample]
        total = sample_variants['pos'].count()
        contamination = sample_variants[sample_variants['filter_status'] == 'contamination']['pos'].count()
        percentage = (contamination / total * 100) if total > 0 else 0
        data.append({'Sample': sample, 'Contamination (%)': percentage, 'Group': group_name})
    return data

ffpe_data = get_contamination_percentages(df_expanded, ffpe_samples, 'FFPE')
smear_data = get_contamination_percentages(df_expanded, smear_samples, 'Smear')

plot_df = pd.DataFrame(ffpe_data + smear_data)

ffpe_contamination = plot_df[plot_df['Group'] == 'FFPE']['Contamination (%)'].tolist()
smear_contamination = plot_df[plot_df['Group'] == 'Smear']['Contamination (%)'].tolist()
mean_ffpe = sum(ffpe_contamination) / len(ffpe_contamination)
mean_smear = sum(smear_contamination) / len(smear_contamination)

t_stat, p_value = ttest_ind(smear_contamination, ffpe_contamination, equal_var=False)

print(f"Mean contamination % in smears: {round(sum(smear_contamination)/len(smear_contamination), 2)}")
print(f"Mean contamination % in FFPE: {round(sum(ffpe_contamination)/len(ffpe_contamination), 2)}")
print(f"t-statistic: {t_stat:.3f}")
print(f"p-value: {p_value:.4f}")
if p_value < 0.05:
    print("Result: The contamination filter is statistically more probable in smears than in FFPE blocks (p < 0.05).")
else:
    print("Result: No statistically significant difference in contamination filter between smears and FFPE blocks (p >= 0.05).")

plt.figure(figsize=(8, 6))
sns.boxplot(x='Group', y='Contamination (%)', data=plot_df)
sns.stripplot(x='Group', y='Contamination (%)', data=plot_df, color='black', jitter=0.2, alpha=0.7, size=6)

plt.text(
    0, max(ffpe_contamination) + 2, f"Mean = {mean_ffpe:.2f}%", 
    ha='center', fontsize=12, color='blue', fontweight='bold'
)
plt.text(
    1, max(smear_contamination) + 2, f"Mean = {mean_smear:.2f}%", 
    ha='center', fontsize=12, color='blue', fontweight='bold'
)

plt.text(
    0.5, max(max(ffpe_contamination), max(smear_contamination)) + 7, 
    f"t = {t_stat:.2f}\np = {p_value:.4f}", 
    ha='center', fontsize=13, color='black', fontweight='bold'
)

plt.title('Contamination Filter Percentage by Sample Type (Boxplot)')
plt.ylabel('Contamination Filter (%)')
plt.xlabel('')
plt.ylim(-15, max(max(ffpe_contamination), max(smear_contamination)) + 15)
plt.tight_layout()
plt.savefig("/mnt/scratch/nazilak/Results/contamination_boxplot.png", dpi=150)
plt.show()


print("Plots saved to /mnt/scratch/nazilak/Results/contamination_boxplot.png")


########################### FIGURE 10.B ##################################

def get_orientation_percentages(df_expanded, sample_set, group_name):
    data = []
    for sample in sample_set:
        sample_variants = df_expanded[df_expanded['sample_nr'].astype(str) == sample]
        total = sample_variants['pos'].count()
        orientation = sample_variants[sample_variants['filter_status'] == 'orientation']['pos'].count()
        percentage = (orientation / total * 100) if total > 0 else 0
        data.append({'Sample': sample, 'orientation (%)': percentage, 'Group': group_name})
    return data

ffpe_data = get_orientation_percentages(df_expanded, ffpe_samples, 'FFPE')
smear_data = get_orientation_percentages(df_expanded, smear_samples, 'Smear')

plot_df = pd.DataFrame(ffpe_data + smear_data)

ffpe_orientation = plot_df[plot_df['Group'] == 'FFPE']['orientation (%)'].tolist()
smear_orientation = plot_df[plot_df['Group'] == 'Smear']['orientation (%)'].tolist()
mean_ffpe = sum(ffpe_orientation) / len(ffpe_orientation)
mean_smear = sum(smear_orientation) / len(smear_orientation)


t_stat, p_value = ttest_ind(smear_orientation, ffpe_orientation, equal_var=False)

print(f"Mean orientation % in smears: {round(sum(smear_orientation)/len(smear_orientation), 2)}")
print(f"Mean orientation % in FFPE: {round(sum(ffpe_orientation)/len(ffpe_orientation), 2)}")
print(f"t-statistic: {t_stat:.3f}")
print(f"p-value: {p_value:.4f}")
if p_value < 0.05:
    print("Result: The orientation filter is statistically more probable in smears than in FFPE blocks (p < 0.05).")
else:
    print("Result: No statistically significant difference in orientation filter between smears and FFPE blocks (p >= 0.05).")


plt.figure(figsize=(8, 6))
sns.boxplot(x='Group', y='orientation (%)', data=plot_df)
sns.stripplot(x='Group', y='orientation (%)', data=plot_df, color='black', jitter=0.2, alpha=0.7, size=6)


plt.text(
    0, max(ffpe_orientation) + 2, f"Mean = {mean_ffpe:.2f}%", 
    ha='center', fontsize=12, color='blue', fontweight='bold'
)
plt.text(
    1, max(smear_orientation) + 2, f"Mean = {mean_smear:.2f}%", 
    ha='center', fontsize=12, color='blue', fontweight='bold'
)

plt.text(
    0.5, max(max(ffpe_orientation), max(smear_orientation)) + 7, 
    f"t = {t_stat:.2f}\np = {p_value:.4f}", 
    ha='center', fontsize=13, color='black', fontweight='bold'
)

plt.title('orientation Filter Percentage by Sample Type (Boxplot)')
plt.ylabel('orientation Filter (%)')
plt.xlabel('')
plt.ylim(-15, max(max(ffpe_orientation), max(smear_orientation)) + 15)
plt.tight_layout()
plt.savefig("/mnt/scratch/nazilak/Results/orientation_boxplot.png", dpi=150)
plt.show()

print("Plots saved to /mnt/scratch/nazilak/Results/ORIENTATION_boxplot.png")