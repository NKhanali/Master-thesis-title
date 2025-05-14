# Comparative Analysis of NGS Tools for Short Somatic Variant Calling
## CLC Genomics Workbench vs Genome Analysis Toolkit

### Abstract

Background: Lung cancer remains one of the leading causes of cancer-related mortality worldwide, with small cell lung cancer (SCLC) accounting for approximately 15% of all cases. Despite advances in diagnostics and treatment, survival rates for SCLC remain poor, underscoring the urgent need for improved therapeutic strategies. Leveraging bioinformatics and DNA sequencing can enhance our understanding of the genetic landscape of SCLCGW and facilitate the development of more effective targeted therapies.

Purpose & aim: In a study investigating the effects of hyperfractionated thoracic radiotherapy on limited-disease SCLCGW patients (THORA), somatic variant calling was initially performed using CLCGW Genomics Workbench (CLCGW), a commercial and proprietary platform known for its user-friendly interface. While CLCGW lowers the barrier for clinical and academic users, its closed-source nature and high cost limit compliance with FAIR (Findable, Accessible, Interoperable, Reusable) data principles. This raised the question of whether an open-source alternative could deliver comparable results. 

Method: To explore this, a somatic variant calling pipeline based on the Genome Analysis Toolkit (GATK) Best Practices was developed and implemented using Snakemake as a workflow manager. This pipeline was applied to the same 52 THORA samples analyzed with CLCGWGW.

Results: Comparative analysis revealed that while CLCGWGW generated a larger set of variant calls, this came at the cost of specificity, as evidenced by a higher proportion of likely false positives. In contrast, GATK's more stringent filtering criteria, though contributing to a lower variant count and reduced concordance with CLCGW, were supported by the biological characteristics of the sample

Conclusion: The GATK pipeline, particularly when combined with workflow managers like Snakemake, offers a scalable, flexible, and fully open-source solution. Although its interface may be less intuitive than that of CLCGW, it compensates with extensive documentation and customization potential, making it a strong long-term option for accurate and reproducible variant analysis.

### Repository Structure
The repository is organized as follows:

- **`Result_figures/`** - The script to all figures produced for the result section. The name of the file corresponds to the figure 
- **`Supplementary_figures/`** - The script to all the supplementary figures. The name of the file corresponds to the figure
- **`CreateHomemadePON.smk/`** - Snakefile scripts for creating the custom panel of normal (PON)
- **`Funcotation_on_CLCdata.smk/`** - Snakefile for running GATK's Funcotation tool on the CLCGW's vcf files
- **`Reference_indexing.smk/`** - Snakefile script for indexing the reference genome
- **`variant_comparison.py/`** - Python script for making a csv file containing information about all the variants called by both CLCGW and GATK
- **`variant_comparison.py/`** - Python script for making a csv file containing information about all the variants called by both CLCGW and GATK, when the custom PON were utilized
- **`Workflow_pipeline_config.yaml/`** - Configuration file for the GATK Best Practices workflow for Short Somatic Variant Discovery
- **'Workflow_pipeline.smk'** - GATK Best Practices workflow for Short Somatic Variant Discovery pipeline using Snakemake workflow manager
