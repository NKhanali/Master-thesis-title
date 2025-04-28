configfile: "/home/nazilak/Master/Script-gatk/COMPLETE_config.yaml"

use_singularity = config.get("use_singularity", True)
singularity_bind = config.get("singularity_bind", "/mnt/scratch:/mnt/scratch")
singularity_image = config.get("singularity_image", "/home/nazilak/gatk_latest.sif")




import glob
import os

# Dynamically detect sample IDs based on input filenames
input_files = glob.glob("/mnt/scratch/nazilak/CLC_vcf/Annotated_somatic_variants/*_B*.vcf")
SAMPLES = [os.path.basename(f).split("_")[0] for f in input_files]

# Rule to define the final output files for all samples (targets funcotator outputs)
rule all:
    input:
        expand("/mnt/scratch/nazilak/CLC_vcf/Funcotator/{sample}_annotatedCLC.vcf", sample=SAMPLES)

# Rule to transform VCF files
rule transform_vcf:
    input:
        lambda wildcards: glob.glob(f"/mnt/scratch/nazilak/CLC_vcf/Annotated_somatic_variants/{wildcards.sample}_B*.vcf")[0]
    output:
        "/mnt/scratch/nazilak/CLC_vcf/Modified/{sample}_modified.vcf"
    shell:
        """
        awk 'BEGIN {{ OFS="\t" }} {{
            if ($0 ~ /^##contig=<ID=MT/) next
            if ($0 ~ /^##contig=<ID=[0-9]+/) {{
                sub(/ID=[0-9]+/, "ID=chr" substr($0, index($0, "ID=") + 3, length($0) - index($0, "ID=") - 1), $0)
            }}
            if ($0 ~ /^##contig=<ID=X/) sub(/ID=X/, "ID=chrX")
            if ($0 ~ /^##contig=<ID=Y/) sub(/ID=Y/, "ID=chrY")

            # Skip ALT == "." variants (but only for data lines, not headers)
            if ($0 !~ /^#/ && $5 == ".") next

            # Prefix chromosome names
            if ($1 ~ /^[0-9]+$/ || $1 == "X" || $1 == "Y") $1 = "chr" $1

            print $0
        }}' {input} > {output}
        """

# Rule to annotate VCF files using Funcotator
rule funcotator:
    input:
        ref="/mnt/scratch/nazilak/resource_bundle_gatk/hg19.fa",
        variants="/mnt/scratch/nazilak/CLC_vcf/Modified/{sample}_modified.vcf",
        data_sources="/mnt/scratch/nazilak/resource_bundle_gatk/data_source"
    output:
        annotated="/mnt/scratch/nazilak/CLC_vcf/Funcotator/{sample}_annotatedCLC.vcf"
    params:
        singularity_cmd=(
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    resources:
        mem_mb=3300
    threads: 1
    shell:
        """
        {params.singularity_cmd} gatk Funcotator \
            -R {input.ref} \
            -V {input.variants} \
            -O {output.annotated} \
            --data-sources-path {input.data_sources} \
            --output-file-format VCF \
            --ref-version hg19 
        """
