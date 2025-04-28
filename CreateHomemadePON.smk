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

use_singularity = config.get("use_singularity", True)
singularity_bind = config.get("singularity_bind", "/mnt/scratch:/mnt/scratch")
singularity_image = config.get("singularity_image", "/home/nazilak/gatk_latest.sif")

use_homemade_pon = config["use_homemade_pon"]



import glob
import os

normals = [
    os.path.splitext(os.path.basename(file))[0].replace("_N_preprocessed", "")
    for file in glob.glob(
        os.path.join(Data_pre_processing, "ApplyBQSR", "*_N_preprocessed.bam")
    )
]

rule mutect2_Normals:
    input:
        ref = Reference_genome + ".fa",
        normal = f"{Data_pre_processing}/ApplyBQSR/{{normal}}_N_preprocessed.bam",
        intervals = f"{Intervals}"
    output:
        normal_pon = f"{Variant_calling}/CreatePON/{{normal}}_pon.vcf.gz"
    params:
        max_mnp_distance = 0,
        Temp_dir = f"{Temp_dir_path}/tmp_M2PON_{{normal}}",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {singularity_bind} {singularity_image} "
            if use_singularity else ""
        )
    log:
        log = f"{logpath}/Mutect2_Normals/{{normal}}.log"
    benchmark: 
        f"{Benchmark}/{{normal}}_Mutect2_PON.txt"
    resources:
        mem_mb = 1500
    threads: 1 
    shell:
        """
        mkdir -p $(dirname {output.normal_pon})
        mkdir -p {params.benchmark_dir} 
        mkdir -p $(dirname {log.log})

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk --java-options "-Xmx4G" Mutect2 \
            -R {input.ref} \
            -I {input.normal} \
            -max-mnp-distance {params.max_mnp_distance} \
            -O {output.normal_pon} \
            2> {log.log}
        """

rule genomicsdb_import:
    input:
        ref = Reference_genome + ".fa",
        pon_vcfs = expand(f"{Variant_calling}/CreatePON/{{normal}}_pon.vcf.gz", normal=normals)
    output:
        pon_db = directory(f"{Variant_calling}/CreatePON/pon_db")  # Flagged as directory
    params:
        intervals = f"{Intervals}",
        Temp_dir = f"{Temp_dir_path}/tmp_GDBI",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {singularity_bind} {singularity_image} "
            if use_singularity else ""
        ),
        pon_vcfs_str = " ".join([f"-V {vcf}" for vcf in expand(f"{Variant_calling}/CreatePON/{{normal}}_pon.vcf.gz", normal=normals)]),
        genomicsdb_cmd = (
             "--update-workspace-path" if config["update_genomicsdb"] else "--genomicsdb-workspace-path"
        )
    log:
        log = f"{logpath}/PON/genomicsdbimport.log"
    benchmark:
            f"{Benchmark}/GenomicsDBImport.txt"
    shell:
        """
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk --java-options '-Xmx8g -Xms8g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport \
            -L {params.intervals} \
            --merge-input-intervals \
            --interval-merging-rule ALL \
            -R {input.ref} \
            {params.genomicsdb_cmd} {output.pon_db} \
            --tmp-dir {params.Temp_dir} \
            {params.pon_vcfs_str} \
            2> {log.log}
        """

rule create_pon:
    input:
        ref = Reference_genome + ".fa",
        db_path = f"{Variant_calling}/CreatePON/pon_db"
    output:
        pon_vcf = f"{Variant_calling}/CreatePON/PanelOfNormals.vcf.gz"
    params:
        Temp_dir = f"{Temp_dir_path}/tmp_CSPON",
        benchmark_dir = f"{Benchmark}",
        singularity_cmd = (
            f"singularity exec --bind {config['singularity']['bind_directory']}:{config['singularity']['bind_directory']} {config['singularity']['singularity_image']} "
            if config["singularity"]["enabled"] else ""
        )    
    log: 
        log = f"{logpath}/PON/create_somatic_pon.log"
    benchmark:
        f"{Benchmark}/CreatePanelOfNormal.txt"
    shell:
        """
        mkdir -p $(dirname {output.pon_vcf})
        mkdir -p {params.benchmark_dir}
        mkdir -p $(dirname {log.log})

        mkdir -p {params.Temp_dir}
        trap "rm -rf {params.Temp_dir}" EXIT

        {params.singularity_cmd} gatk \
            CreateSomaticPanelOfNormals \
            -R {input.ref} \
            -V gendb://{input.db_path} \
            -O {output.pon_vcf} \
            2> {log.log}
        """

    