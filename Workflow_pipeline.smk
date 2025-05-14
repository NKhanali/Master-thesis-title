configfile: "/home/nazilak/Master/Script-gatk/COMPLETE_config.yaml"

Reference_genome = config["Reference_genome"]
logpath = config["logpath"]
Data_pre_processing = config["Data_pre_processing"]
Variant_calling = config["Variant_calling"]
Temp_dir_path = config["Temp_dir_path"]
Resource_bundle_gatk = config["Resource_bundle_gatk"]
homemade_PON = config["homemade_PON"]
Script = config["Script"]
Intervals = config["Intervals"]
Benchmark = config["Benchmark"]
sample_sheet = config["sample_sheet"] 
CreateHomemadePON = config["CreateHomemadePON"]
use_singularity = config.get("use_singularity", True)
singularity_bind = config.get("singularity_bind", "/mnt/scratch:/mnt/scratch")
singularity_image = config.get("singularity_image", "/home/nazilak/gatk_latest.sif")

use_homemade_pon = config["use_homemade_pon"]


import pandas as pd 
read_sample_sheet = pd.read_csv(sample_sheet, dtype="str").set_index("Sample", drop=False)

def NR1_sample_sheet(wildcards):
    return read_sample_sheet.loc[wildcards.Sample, "Normal_R1"]

def NR2_sample_sheet(wildcards):
    return read_sample_sheet.loc[wildcards.Sample, "Normal_R2"]

def TR1_sample_sheet(wildcards):
    return read_sample_sheet.loc[wildcards.Sample, "Tumor_R1"]

def TR2_sample_sheet(wildcards):
    return read_sample_sheet.loc[wildcards.Sample, "Tumor_R2"]

import csv 
def get_sample_number(csv_file):
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        samples = [row["Sample"] for row in reader if row["Sample"]]
    return samples

SAMPLES_NUMBER = get_sample_number(sample_sheet)


if use_homemade_pon:
    include: f"{CreateHomemadePON}/CreateHomemadePON.smk"

rule all:
    input: 
        expand(
            f"{homemade_PON}/Liftover_PON/{{Sample}}_liftover.vcf", Sample=SAMPLES_NUMBER
        ) if use_homemade_pon else expand(
            f"{Variant_calling}/Liftover/{{Sample}}_liftover.vcf", Sample=SAMPLES_NUMBER
        ),

    

rule Bwa_alignment_normal:
    input:   
        NR1 = NR1_sample_sheet,
        NR2 = NR2_sample_sheet,
        ref = Reference_genome + ".fa"  
    output: 
        Sam = f"{Data_pre_processing}/Alignment_normal/{{Sample}}_N.sam"
    params:
        benchmark_dir = f"{Benchmark}"
    log:    
        log = f"{logpath}/Alignment/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_bwa_alignment_normal.txt"
    resources:
        mem_mb = 6450
    shell:
        """
        mkdir -p $(dirname {output.Sam})
        mkdir -p {params.benchmark_dir}  
        mkdir -p $(dirname {log.log})

        bwa mem -t 10 -R '@RG\\tID:{wildcards.Sample}_N\\tSM:{wildcards.Sample}_N\\tPL:illumina\\tLB:normal' {input.ref} {input.NR1} {input.NR2} > {output.Sam} 2> {log.log}
        """

rule Sam_to_bam_normal:
    input:
        Sam= f"{Data_pre_processing}/Alignment_normal/{{Sample}}_N.sam",
        interval = f"{Intervals}"
    output:
        Bam= f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N.bam"
    params:
        benchmark_dir = f"{Benchmark}"
    benchmark: 
        f"{Benchmark}/{{Sample}}_sam_to_bam_normal.txt"
    resources:
        mem_mb = 6
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.Bam})
        mkdir -p {params.benchmark_dir}

        samtools view -Sb {input.Sam} -L {input.interval} -o {output.Bam} 
        """
rule Sort_normal:
    input: 
        Bam = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N.bam"
    output: 
        sort = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N_sorted.bam",
        index = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N_sorted.bam.bai"
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_sortN_{{Sample}}",
        benchmark_dir = f"{Benchmark}"
    log:
        log = f"{logpath}/Sam_to_bam/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_sort_normal.txt"
    resources:
        mem_mb = 900  
    threads: 1 
    shell: 
        """
        mkdir -p $(dirname {output.sort})
        mkdir -p $(dirname {log.log}) 
        mkdir -p {params.benchmark_dir}

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        samtools sort -n -@ {threads} -m {resources.mem_mb}M  -T {params.Temp_dir}/tmp_prefix {input.Bam} -o {output.sort} 2> {log.log}
        samtools index {output.sort}
        """
        

rule MarkDuplicateSpark_Normal:
    input: 
        normal_bam = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N_sorted.bam",
        intervals = f"{Intervals}"
    output: 
        marked_normal_bam = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_N_marked.bam"
    params:
        metrics_normal = f"{Data_pre_processing}/MarkDuplicateSpark/metrics/{{Sample}}_N_metrics.txt",
        Temp_dir = f"{Temp_dir_path}/tmp_MD_{{Sample}}N",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/MarkDuplicateSpark/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_MarkDuplicateSpark_Normal.txt"
    resources:
        mem_mb = 9850
    shell:
        """
        mkdir -p $(dirname {output.marked_normal_bam})
        mkdir -p $(dirname {params.metrics_normal})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd}gatk --java-options "-Xmx8G -Djava.io.tmpdir={params.Temp_dir}" MarkDuplicatesSpark \
            --conf 'spark.local.dir={params.Temp_dir}' \
            --conf 'spark.executor.memory=8g' \
            --conf 'spark.driver.memory=8g' \
            --conf 'spark.executor.cores=20' \
            -I {input.normal_bam} \
            -L {input.intervals} \
            -O {output.marked_normal_bam} \
            --tmp-dir {params.Temp_dir} \
            -M {params.metrics_normal} \
            2> {log.log} \
            --spark-verbosity WARN \
            --verbosity WARNING
        """

rule BaseRecalibrator_Normal:
    input: 
        ref = Reference_genome + ".fa",
        MNB = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_N_marked.bam",
        database = f"{Resource_bundle_gatk}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
        database_index = f"{Resource_bundle_gatk}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        intervals = f"{Intervals}"
    output: 
        Recalibrated_table_normal = f"{Data_pre_processing}/BaseRecalibrator/{{Sample}}_N_recal_data.table"
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_BR_{{Sample}}N",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/BaseRecalibrator/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_BaseRecalibrator_Normal.txt"   
    resources:
        mem_mb = 1650
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.Recalibrated_table_normal})
        mkdir -p {params.Temp_dir}
        mkdir -p $(dirname {log.log})  
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk BaseRecalibrator \
            -L {input.intervals} \
            -I {input.MNB} \
            -R {input.ref} \
            --known-sites {input.database} \
            -O {output.Recalibrated_table_normal} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """

rule ApplyBQSR_normal:
    input: 
        ref = Reference_genome + ".fa", 
        intervals = f"{Intervals}",
        MNB = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_N_marked.bam",
        NRT = f"{Data_pre_processing}/BaseRecalibrator/{{Sample}}_N_recal_data.table"
    output: 
        preprocessed_normal = f"{Data_pre_processing}/ApplyBQSR/{{Sample}}_N_preprocessed.bam"
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_ABQSR_{{Sample}}N",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/ApplyBQSR/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_ApplyBQSR_normal.txt"
    resources:
        mem_mb = 1250
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.preprocessed_normal})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})  
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk ApplyBQSR \
            -R {input.ref} \
            -L {input.intervals} \
            -I {input.MNB} \
            --bqsr-recal-file {input.NRT} \
            -O {output.preprocessed_normal} \
            2> {log.log} \
            --tmp-dir {params.Temp_dir} \
        """
rule GetPileupSummaries_normal:
    input: 
        normal = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N_sorted.bam",
        normal_index = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_N_sorted.bam.bai",
        germline_resource = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz",
        germline_resource_index = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz.tbi",
        intervals = f"{Intervals}"
    output: 
        summaries_normal = f"{Variant_calling}/GetPileupSummaries/{{Sample}}_N.table"
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_GPS_{{Sample}}N",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/GetPileupSummaries/{{Sample}}N.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_GetPileupsummaries_normal.txt"
    threads: 1 
    resources:
        mem_mb = 2000
    shell:
        """
        mkdir -p $(dirname {output.summaries_normal})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk  GetPileupSummaries\
            -I {input.normal} \
            -V {input.germline_resource} \
            -L {input.intervals} \
            -O {output.summaries_normal} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """


rule Bwa_alignment_tumor:
    input:
        TR1 = TR1_sample_sheet,
        TR2 = TR2_sample_sheet,
        ref = Reference_genome + ".fa"               
    output:
        Sam = f"{Data_pre_processing}/Alignment_tumor/{{Sample}}_T.sam"
    params:
        benchmark_dir = f"{Benchmark}"
    log:
        log = f"{logpath}/Alignment/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_bwa_alignment_tumor.txt"
    resources:
        mem_mb = 6450
    shell:
        """
        mkdir -p $(dirname {output.Sam})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})  
        bwa mem -t 10 -R '@RG\\tID:{wildcards.Sample}_T\\tSM:{wildcards.Sample}_T\\tPL:illumina\\tLB:normal' {input.ref} {input.TR1} {input.TR2} > {output.Sam} 2> {log.log}
        """

rule Sam_to_bam_tumor:
    input:
        Sam= f"{Data_pre_processing}/Alignment_tumor/{{Sample}}_T.sam",
        intervals = f"{Intervals}"
    output:
        Bam = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T.bam"
    params: 
        benchmark_dir = f"{Benchmark}"
    benchmark: 
        f"{Benchmark}/{{Sample}}_sam_to_bam_tumor.txt"
    resources:
        mem_mb= 6
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.Bam})
        mkdir -p {params.benchmark_dir}
    
        samtools view -Sb {input.Sam} -L {input.intervals} -o {output.Bam}
        """

rule Sort_tumor:
    input:
        Bam = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T.bam" 
    output: 
        sort = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T_sorted.bam",
        index = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T_sorted.bam.bai"
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_sortT_{{Sample}}",
        benchmark_dir = f"{Benchmark}"
    log:
        log = f"{logpath}/Sam_to_bam/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_sort_tumor.txt"
    resources:
        mem_mb = 900
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.sort})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        samtools sort -@ {threads} -m {resources.mem_mb}M -T {params.Temp_dir}/tmp_prefix {input.Bam} -o {output.sort} 2> {log.log}
        samtools index {output.sort}
        """ 
    

rule MarkDuplicateSpark_Tumor:
    input:
        tumor_bam = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T_sorted.bam",
        intervals = f"{Intervals}"
    output:
        marked_tumor_bam = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_T_marked.bam"
    params:
        metrics_tumor = f"{Data_pre_processing}/MarkDuplicateSpark/metrics/{{Sample}}_T_metrics.txt",
        benchmark_dir = f"{Benchmark}",
        Temp_dir = f"{Temp_dir_path}/tmp_MD_{{Sample}}T",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/MarkDuplicateSpark/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_MarkDuplicateSpark_Tumor.txt"
    resources:
        mem_mb = 9850 
    shell:
        """
        mkdir -p $(dirname {output.marked_tumor_bam})   
        mkdir -p $(dirname {params.metrics_tumor})
        mkdir -p $(dirname {log.log}) 
        mkdir -p {params.benchmark_dir}
         
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk --java-options "-Xmx8G -Djava.io.tmpdir={params.Temp_dir}" MarkDuplicatesSpark \
            --conf 'spark.local.dir={params.Temp_dir}' \
            --conf 'spark.executor.memory=8g' \
            --conf 'spark.driver.memory=8g' \
            --conf 'spark.executor.cores=20' \
            -I {input.tumor_bam} \
            -L {input.intervals} \
            -O {output.marked_tumor_bam} \
            --tmp-dir {params.Temp_dir} \
            -M {params.metrics_tumor} \
            2> {log.log} \
            --spark-verbosity WARN \
            --verbosity WARNING 
        """

rule BaseRecalibrator_Tumor:
    input: 
        ref = Reference_genome + ".fa",
        MTB = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_T_marked.bam",
        intervals = f"{Intervals}",
        database = f"{Resource_bundle_gatk}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
        database_index = f"{Resource_bundle_gatk}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.idx"
    output: 
        Recalibrated_table_tumor = f"{Data_pre_processing}/BaseRecalibrator/{{Sample}}_T_recal_data.table"
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_BR_{{Sample}}T",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/BaseRecalibrator/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_BaseRecalibrator_Tumor.txt"
    resources:
        mem_mb = 1650
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.Recalibrated_table_tumor})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})  
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk BaseRecalibrator \
            -L {input.intervals} \
            -I {input.MTB} \
            -R {input.ref} \
            --known-sites {input.database} \
            -O {output.Recalibrated_table_tumor} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """

rule ApplyBQSR_tumor:
    input: 
        ref = Reference_genome + ".fa", 
        intervals = f"{Intervals}",
        MTB = f"{Data_pre_processing}/MarkDuplicateSpark/{{Sample}}_T_marked.bam",
        TRT = f"{Data_pre_processing}/BaseRecalibrator/{{Sample}}_T_recal_data.table"
    output: 
        preprocessed_tumor = f"{Data_pre_processing}/ApplyBQSR/{{Sample}}_T_preprocessed.bam",
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_ABQSR_{{Sample}}T",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/ApplyBQSR/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_ApplyBQSR_tumor.txt"
    resources:
        mem_mb = 1250
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.preprocessed_tumor})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})  

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT
        
        {params.singularity_cmd} gatk ApplyBQSR \
            -R {input.ref} \
            -L {input.intervals} \
            -I {input.MTB} \
            --bqsr-recal-file {input.TRT} \
            -O {output.preprocessed_tumor} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """

rule GetPileupSummaries_tumor:
    input: 
        tumor = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T_sorted.bam",
        tumor_index = f"{Data_pre_processing}/Sam_to_bam/{{Sample}}_T_sorted.bam.bai",
        germline_resource = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz",
        germline_resource_index = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz.tbi",
        intervals = f"{Intervals}"
    output: 
        summaries_tumor = f"{Variant_calling}/GetPileupSummaries/{{Sample}}_T.table"
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_GPS_{{Sample}}T",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/GetPileupSummaries/{{Sample}}T.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_GetPileupSummaries_tumor.txt"
    threads: 1 
    resources:
        mem_mb = 2000
    shell:
        """
        mkdir -p $(dirname {output.summaries_tumor})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT
    
        {params.singularity_cmd} gatk  GetPileupSummaries\
            -I {input.tumor} \
            -V {input.germline_resource} \
            -L {input.intervals} \
            -O {output.summaries_tumor} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """



rule CalculateContamination:
    input: 
        summaries_tumor = f"{Variant_calling}/GetPileupSummaries/{{Sample}}_T.table",
        summaries_normal = f"{Variant_calling}/GetPileupSummaries/{{Sample}}_N.table"
    output: 
        contamination = f"{Variant_calling}/CalculateContamination/{{Sample}}_contamination.table",
        tumor_segmentation = f"{Variant_calling}/CalculateContamination/{{Sample}}_tumor_segments.tsv"
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_CC_{{Sample}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = f"{logpath}/CalculateContamination/{{Sample}}.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_CalculateContamination.txt"
    resources:
        mem_mb = 400
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.contamination})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk CalculateContamination\
            -I {input.summaries_tumor} \
            -matched {input.summaries_normal} \
            -O {output.contamination} \
            --tumor-segmentation {output.tumor_segmentation} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """





rule Mutect2:
    input:
        ref = Reference_genome + ".fa",
        tumor = f"{Data_pre_processing}/ApplyBQSR/{{Sample}}_T_preprocessed.bam",
        normal = f"{Data_pre_processing}/ApplyBQSR/{{Sample}}_N_preprocessed.bam",
        intervals = f"{Intervals}",
        germline_resource = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz",
        germline_resource_index = f"{Resource_bundle_gatk}/af-only-gnomad.hg38.vcf.gz.tbi",
        pon =
            f"{Variant_calling}/CreatePON/PanelOfNormals.vcf.gz"
            if use_homemade_pon else 
            f"{Resource_bundle_gatk}/1000g_pon.hg38.vcf.gz",
    output:
        somatic_vcf = 
            f"{homemade_PON}/Mutect2_PON/{{Sample}}.vcf.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/Mutect2/{{Sample}}.vcf.gz",
        f1r2_tar_gz =
            f"{homemade_PON}/Mutect2_PON/{{Sample}}.f1r2.tar.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/Mutect2/{{Sample}}.f1r2.tar.gz",     
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_M2_{{Sample}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log =
            f"{logpath}/Mutect2_PON/{{Sample}}.log"
            if use_homemade_pon else 
            f"{logpath}/Mutect2/{{Sample}}.log",
    benchmark:
        f"{Benchmark}/{{Sample}}_Mutect2.txt" if not use_homemade_pon else f"{Benchmark}/{{Sample}}_Mutect2_PON.txt"
    resources:
        mem_mb = 1500
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.somatic_vcf})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk --java-options "-Xmx4G" Mutect2 \
            -R {input.ref} \
            -L {input.intervals} \
            -I {input.tumor} \
            -I {input.normal} \
            -normal {wildcards.Sample}_N \
            --germline-resource {input.germline_resource} \
            --panel-of-normals {input.pon} \
            -O {output.somatic_vcf} \
            --f1r2-tar-gz {output.f1r2_tar_gz} \
            --af-of-alleles-not-in-resource 1e-6 \
            --tmp-dir {params.Temp_dir} \
            2> {log.log} \
        """

rule LearnReadOrientationModel:
    input: 
        f1r2_tar_gz =
            
            f"{homemade_PON}/Mutect2_PON/{{Sample}}.f1r2.tar.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/Mutect2/{{Sample}}.f1r2.tar.gz",
    output: 
        artifacts =
            f"{homemade_PON}/LearnReadOrientationModel_PON/{{Sample}}.artifacts.tar.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/LearnReadOrientationModel/{{Sample}}.artifacts.tar.gz",
            
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_LROM_{{Sample}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )

    log:
        log =
            f"{logpath}/LearnReadOrientationModel_PON/{{Sample}}.log"
            if use_homemade_pon else 
            f"{logpath}/LearnReadOrientationModel/{{Sample}}.log",
    benchmark:  
        benchmark =
            f"{Benchmark}/{{Sample}}_LearnReadOrientationModel_PON.txt"
            if use_homemade_pon else 
            f"{Benchmark}/{{Sample}}_LearnReadOrientationModel.txt",
    resources:
        mem_mb = 1000
    threads: 1

    shell:
        """
        mkdir -p $(dirname {output.artifacts})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk --java-options "-Xmx4G" LearnReadOrientationModel \
            -I {input.f1r2_tar_gz} \
            -O {output.artifacts} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """
 

rule FilterMutectCalls:
    input: 
        ref = Reference_genome + ".fa",
        vcf = 
            f"{homemade_PON}/Mutect2_PON/{{Sample}}.vcf.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/Mutect2/{{Sample}}.vcf.gz",
            
        contamination_table = f"{Variant_calling}/CalculateContamination/{{Sample}}_contamination.table",
        tumor_segmentation = f"{Variant_calling}/CalculateContamination/{{Sample}}_tumor_segments.tsv",
        artifacts = 
            f"{homemade_PON}/LearnReadOrientationModel_PON/{{Sample}}.artifacts.tar.gz"
            if use_homemade_pon else 
            f"{Variant_calling}/LearnReadOrientationModel/{{Sample}}.artifacts.tar.gz",
        intervals = f"{Intervals}",
    output: 
        filtered_variants = 
            f"{homemade_PON}/FilterMutectCalls_PON/{{Sample}}_filtered_variants.vcf.gz"
            if use_homemade_pon else
            f"{Variant_calling}/FilterMutectCalls/{{Sample}}_filtered_variants.vcf.gz"
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_FMC_{{Sample}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = 
            f"{logpath}/FilterMutectCalls_PON/{{Sample}}.log"
            if use_homemade_pon else
            f"{logpath}/FilterMutectCalls/{{Sample}}.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_FilterMutectCalls_PON.txt"
        if use_homemade_pon else
        f"{Benchmark}/{{Sample}}_FilterMutectCalls.txt"
        
    resources:
        mem_mb = 1200
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.filtered_variants})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})  
        
        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk FilterMutectCalls \
            -R {input.ref} \
            -L {input.intervals} \
            -V {input.vcf} \
            --contamination-table {input.contamination_table} \
            --tumor-segmentation {input.tumor_segmentation} \
            --orientation-bias-artifact-priors {input.artifacts} \
            -O {output.filtered_variants} \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """

rule Funcotator:
    input: 
        ref = Reference_genome + ".fa",
        variants = 
            f"{homemade_PON}/FilterMutectCalls_PON/{{Sample}}_filtered_variants.vcf.gz"
            if use_homemade_pon else
            f"{Variant_calling}/FilterMutectCalls/{{Sample}}_filtered_variants.vcf.gz",
        data_sources = f"{Resource_bundle_gatk}/data_source",
        intervals = f"{Intervals}"
    output: 
        annotated = 
            f"{homemade_PON}/Funcotator_PON/{{Sample}}_annotated.vcf"
            if use_homemade_pon else
            f"{Variant_calling}/Funcotator/{{Sample}}_annotated.vcf"
 
    params: 
        Temp_dir = f"{Temp_dir_path}/tmp_FCT_{{Sample}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = 
            f"{logpath}/Funcotator_PON/{{Sample}}.log"
            if use_homemade_pon else
            f"{logpath}/Funcotator/{{Sample}}.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_Funcotator_PON.txt"
        if use_homemade_pon else
        f"{Benchmark}/{{Sample}}_Funcotator.txt"
    resources:
        mem_mb = 3300
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.annotated})
        mkdir -p $(dirname {log.log})  
        mkdir -p {params.benchmark_dir}      

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk Funcotator \
            -R {input.ref} \
            -V {input.variants} \
            -L {input.intervals} \
            -O {output.annotated} \
            --data-sources-path {input.data_sources} \
            --output-file-format VCF \
            --ref-version hg38 \
            --tmp-dir {params.Temp_dir} \
            2> {log.log}
        """

rule Liftover:
    input: 
        refHG19 = f"{Resource_bundle_gatk}/hg19.fa",
        chain38to19 = f"{Resource_bundle_gatk}/hg38ToHg19.over.chain.gz",
        vcfGATK = 
            f"{homemade_PON}/Funcotator_PON/{{Sample}}_annotated.vcf"
            if use_homemade_pon else
            f"{Variant_calling}/Funcotator/{{Sample}}_annotated.vcf"
    output: 
        lifted =  
            f"{homemade_PON}/Liftover_PON/{{Sample}}_liftover.vcf"
            if use_homemade_pon else
            f"{Variant_calling}/Liftover/{{Sample}}_liftover.vcf",
        Reject = 
            f"{homemade_PON}/Liftover_PON/{{Sample}}_reject.vcf"
            if use_homemade_pon else
            f"{Variant_calling}/Liftover/{{Sample}}_reject.vcf",
    params: 
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )
    log:
        log = 
            f"{logpath}/Liftover_PON/{{Sample}}.log"
            if use_homemade_pon else
            f"{logpath}/Liftover/{{Sample}}.log"
    benchmark: 
        f"{Benchmark}/{{Sample}}_Liftover_PON.txt"
        if use_homemade_pon else
        f"{Benchmark}/{{Sample}}_Liftover.txt"
    resources:
        mem_mb = 3800
    threads: 1
    shell: 
        """
        mkdir -p $(dirname {output.lifted})
        mkdir -p {params.benchmark_dir}   
        mkdir -p $(dirname {log.log})
               
        {params.singularity_cmd} gatk LiftoverVcf \
            -I {input.vcfGATK} \
            -R {input.refHG19} \
            -C {input.chain38to19} \
            -O {output.lifted} \
            --REJECT {output.Reject} \
            --RECOVER_SWAPPED_REF_ALT \
            2> {log.log}
        """


