import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns


def vcf(file_path, workflow, sample_type=None):
    data = []
    samples = []

    with open(file_path, 'r', encoding='latin1') as f:
        for line in f:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = list(set([sample.split('_')[0] for sample in parts[9:]]))
            elif not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom = parts[0]
                if workflow == 'CLC' and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                info = parts[7]

                if workflow == 'GATK':
                    filt = parts[6]
                elif workflow == 'CLC':
                    filt = sample_type

                
                gene = ''
                variant_class = ''
                if "FUNCOTATION=" in info:
                    try:
                        funcotation_str = info.split("FUNCOTATION=")[1]
                        funcotation_str = funcotation_str.strip('[]')
                        funcotation_entries = funcotation_str.split('],[')
                        for entry in funcotation_entries:
                            entry = entry.strip('[]')
                            funcotation_values = entry.split('|')
                            if len(funcotation_values) > 5:
                                gene_candidate = funcotation_values[0].strip().upper()
                                variant_class_candidate = funcotation_values[5].strip().upper()
                                
                                if gene_candidate and variant_class_candidate:
                                    gene = gene_candidate
                                    variant_class = variant_class_candidate
                                    break  
                    except Exception as e:
                        print(f"Warning: Failed to parse FUNCOTATION in line with error: {e}")

                for sample in samples:
                    data.append([
                        chrom,
                        int(pos),
                        ref,
                        alt.split(',')[0],
                        sample,
                        workflow,
                        filt,
                        gene,
                        variant_class
                    ])


    return pd.DataFrame(data, columns=[
        'chr', 'pos', 'ref', 'alt',
        'sample_nr', 'workflow', 'filter_status',
        'gene', 'variant_class'
    ])




gatk_path = "/mnt/scratch/nazilak/GATK/Variant_calling/Liftover/"
clc_annotated_path = "/mnt/scratch/nazilak/CLC_vcf/Funcotator/"
clc_low_freq_path = "/mnt/scratch/nazilak/CLC_vcf/Somatic_variants_with_low_frequency_in_germline/"


def process_files(gatk_path, clc_annotated_path, clc_low_freq_path):
    all_data = []

    for vcf_file in glob.glob(os.path.join(gatk_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "GATK")
        df['sample_nr'] = sample_name
        all_data.append(df)

    for vcf_file in glob.glob(os.path.join(clc_annotated_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "annotated")
        df['sample_nr'] = sample_name
        all_data.append(df)

    for vcf_file in glob.glob(os.path.join(clc_low_freq_path, "*.vcf")):
        sample_name = os.path.basename(vcf_file).split('_')[0]
        df = vcf(vcf_file, "CLC", "low_frequency")
        df['sample_nr'] = sample_name
        all_data.append(df)


    combined = pd.concat(all_data, ignore_index=True)

    duplicates_within_sample = combined[combined.duplicated(subset=[
        'chr', 'pos', 'ref', 'alt', 'sample_nr', 'workflow', 'filter_status'
    ], keep=False)]
    print("Duplicate rows found within the same sample:")
    print(duplicates_within_sample)

    return combined

combined_data = process_files(gatk_path, clc_annotated_path, clc_low_freq_path)
combined_data.to_csv("/mnt/scratch/nazilak/Results/variant_comparison.csv", index=False)
print("Combined data saved to /mnt/scratch/nazilak/Results/variant_comparison.csv")

