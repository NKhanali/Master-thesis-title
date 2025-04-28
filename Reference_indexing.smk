
Reference_genome = "/mnt/scratch/nazilak/ref_index/reference_hg38_with_chr.fa"

rule all:
    input: 
        expand("/mnt/scratch/nazilak/resource_bundle_gatk/GCF_000001405.13_GRCh37_genomic.{ext}", ext = ["amb", "ann", "bwt", "pac", "sa", "dict", "fna.fai"])


rule ref_index:
    input: 
        ref = Reference_genome + ".fna"
    output:
        Reference_genome + ".amb",  
        Reference_genome + ".ann",  
        Reference_genome + ".bwt",  
        Reference_genome + ".pac",  
        Reference_genome + ".sa",
        #test med minichromosom 
    shell: 
        """
        bwa index {input.ref} 
        """
        # Add -p option to specify prefix -p {Reference_genome}


rule fai:
    input:
        ref = Reference_genome + ".fna"  
    output: 
        Reference_genome + ".fna.fai"
    shell:
        """
        samtools faidx {input.ref} -o {output}
        echo "Generated {output}" && ls -lh {output}
        """
rule dict:
    input:
        ref = Reference_genome + ".fna"
    output: 
        Reference_genome + ".dict"
    shell:
        """
        singularity exec --bind /mnt/scratch:/mnt/scratch /home/nazilak/gatk_latest.sif \
        gatk CreateSequenceDictionary -R {input.ref} -O {output}
        echo "Generated {output}" && ls -lh {output}
        """