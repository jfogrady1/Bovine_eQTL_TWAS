import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes

ALL_genes_eQTL = pd.read_csv("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", sep = "\t")
print(ALL_genes_eQTL.head())
ALL_genes_gcta = ALL_genes_eQTL["phenotype_id"].tolist()
CONTROL_genes_eQTL = pd.read_csv("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", sep = "\t")
CONTROL_genes_gcta = CONTROL_genes_eQTL["phenotype_id"].tolist()
INFECTED_genes_eQTL = pd.read_csv("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", sep = "\t")
INFECTED_genes_gcta = CONTROL_genes_eQTL["phenotype_id"].tolist()




rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz"),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/selected.{sample}.filtered.vcf.gz", sample = config["samples"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}_rename.txt", sample = config["samples"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_RNA_ALL_imputation_summary_statistics_sorted.info",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/continuous_covar.qcovar",
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/add_h2/{gene}.phen", gene = ALL_genes_gcta),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/Cis-h2-figure.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/plots/Imputation_performance_DNA_V_DNA_RNA.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.FILTERED.dose.vcf.gz",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/private_RNA.ld",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/PNPLA1_ieQTL.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/AD_LRT/ALL_AD_results.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/plots/Gprofiler_ALL_DE_genes_jitter.pdf"



##################################################
# ieQTL mapping
##################################################

rule qtl_mapping_interaction:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi"),# Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt", # Covariates
        interaction_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_TB_interaction_input.txt"
    output:
        outputdir="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor", # Genotypes
       prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction"
    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/
        python3 -m tensorqtl {params.plink_prefix_path} {input.bed_4_tensor[0]} {params.prefix} \
            --covariates {input.covariates_file} \
            --window 1000000 \
            --seed 1864 \
            --interaction {input.interaction_file} \
            --maf_threshold_interaction 0.05 \
            --mode cis_nominal
        '''


# ieQTL plotting
rule ieQTL_plot:
    input:
        ieQTL_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz",
        vcf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz",
        expression = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_residualised_expression.txt',
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/plot_ieQTL.R"
    output:
        GBP4 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/GBP4_ieQTL.pdf",
        PCBP2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/PCBP2_ieQTL.pdf",
        RELT = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/RELT_ieQTL.pdf",
        PLD3 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/PLD3_ieQTL.pdf",
        PNPLA1 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/PNPLA1_ieQTL.pdf",

    shell:
        '''
        Rscript {input.script} {input.ieQTL_results} {input.vcf} {input.expression} {output.GBP4} {output.PCBP2} {output.RELT} {output.PLD3} {output.PNPLA1}
        '''


#########################################################
# Variant calling from RNA-seq data
#########################################################

# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz

# 1. Add read groups, sort, mark duplicates, and create index
rule task1:
    input:
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        read_group = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/rg_added_{sample}-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/dedupped_{sample}-STARAligned.sortedByCoord.out.bam",
        duplicated_metrics = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/{sample}_marked_dup_metrix.txt"
    threads: 30
    params: TMPDIR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/"
    shell:
        '''        
        java -jar {input.picard} AddOrReplaceReadGroups I={input.in_bam} O={output.read_group} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        {input.gatk} MarkDuplicates -I {output.read_group} -O {output.dedupped} \
        -CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M {output.duplicated_metrics}
        '''

# 2. Split N trim and reassign mapping qualities
rule task2:
    input:
        reference="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        dedupped ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/dedupped_{sample}-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/"
    threads: 10
    shell:
        '''
        {input.gatk} SplitNCigarReads -R {input.reference} -I {input.dedupped} -O {output.split} 
        '''


# 3. Base recalibration (BQSR)
rule task3:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        cigar = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam"
    output:
        recalibration_table="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/BQSR/{sample}-recal.table",
        BQSR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/BQSR/BQSR_{sample}-STARAligned.sortedByCoord.out.bam"

    shell:
        '''

        {input.gatk}  BaseRecalibrator -I {input.cigar} -R {input.reference} --known-sites {input.dbSNP} -O {output.recalibration_table}
        {input.gatk}  ApplyBQSR -R {input.reference} -I {input.cigar} --bqsr-recal-file {output.recalibration_table} -O {output.BQSR}
        '''


# 4. Run the haplotypecaller
rule task4:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        BQSR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/BQSR/BQSR_{sample}-STARAligned.sortedByCoord.out.bam",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
    
    output:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}.vcf.gz"

    shell:
        '''
        {input.gatk} HaplotypeCaller -R {input.reference} -I {input.BQSR} -O {output.vcf} -dbsnp {input.dbSNP} --dont-use-soft-clipped-bases --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 0
        '''

rule filter_snps:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",

    output:
        out_filter="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}.filtered.vcf.gz"
    shell:
        " {input.gatk} " +
        " VariantFiltration " +
        " -R {input.reference} " +
        " -V {input.vcf} " +
        " -O {output.out_filter} " +
        " --verbosity ERROR " +
        " --filter \"QD < 2.0\" --filter-name \"snp_QD2\" "
        " --filter \"FS > 30.0\" --filter-name \"snp_FS30\" " 
        " --filter \"DP < 4.0\" --filter-name \"snp_DP4\" "




# 6. Selected variants
rule select_variants:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        in_filter="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}.filtered.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk"

    output:
        out_filter2="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/selected.{sample}.filtered.vcf.gz",
        
    shell:
        '''
        {input.gatk} SelectVariants -R {input.reference} -V {input.in_filter} --select-type-to-include SNP --exclude-filtered -O {output.out_filter2}
        '''

rule rename_vcf:
    input:
        out_filter2="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/selected.{sample}.filtered.vcf.gz"

    output:
        out_filter3="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/renamed/{sample}.selected.vcf.gz",
        name_file ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}_rename.txt"

    params:
        name = "{sample}"
    
    shell:
        """
        echo "20 {params.name}" > {output.name_file}
        bcftools reheader -s {output.name_file} {input.out_filter2} -o {output.out_filter3} && tabix -p vcf {output.out_filter3}
        """

rule merge_RNA_vcf:
    input: 
        vcfs = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/renamed/{sample}.selected.vcf.gz", sample = config["samples"])
    
    output:
        merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_merged.vcf.gz"

    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.merged} && tabix -p vcf {output.merged}
        """

rule filter_RNA_seq_variant_call:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_merged.vcf.gz"

    output:
        vcf_filtered = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_filtered.vcf.gz"

    shell:
        """
        bcftools view -Oz -m2 -M2 -v snps -i 'F_MISSING < 0.2' -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 {input.vcf} -o {output.vcf_filtered} && tabix -p vcf {output.vcf_filtered}
        """

rule change_SNP_data_name:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_Final.vcf.gz",
        names = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Rename.txt",
        header = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/header.txt"

    output:
        vcf_bcftools =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_bcftools_format.vcf.gz",
        vcf_renamed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_Final_renamed_temp.vcf.gz",
        vcf_renamed_header = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_Final_renamed.vcf.gz",
        vcf_renamed_header_noPR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_Final_renamed_noPR.vcf.gz"

    shell:
        """
        # correct bcftools format
        bcftools view  {input.vcf} -Oz -o {output.vcf_bcftools}
        bcftools reheader -s {input.names} {output.vcf_bcftools} -o {output.vcf_renamed}
        bcftools reheader -h {input.header} {output.vcf_bcftools} -o {output.vcf_renamed_header}
        bcftools annotate -x  "INFO/PR" {output.vcf_renamed_header} -Oz -o {output.vcf_renamed_header_noPR} && tabix -p vcf {output.vcf_renamed_header_noPR}
        """

rule concat_with_SNP_data:
    input:
        RNA_seq = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_filtered.vcf.gz",
        SNP_data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_Final_renamed_noPR.vcf.gz"

    output:
        RNA_SNP_concat = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_data_RNA_concat_final.vcf.gz"

    shell:
        """
        bcftools concat -a -d snps {input.SNP_data} {input.RNA_seq} -Oz -o {output.RNA_SNP_concat} && tabix -p vcf {output.RNA_SNP_concat}
        """
    


rule split_target_RNA_SNP_by_chr:
    input:
        cleaned_target = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_data_RNA_concat_final.vcf.gz"
    output:
        target_chromosome = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_CHR{chromosome}.vcf.gz",
    params:
        chromosome = "{chromosome}",
        path = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_CHR",
        ext = ".vcf.gz"
    shell:
        '''
        vcftools --gzvcf {input.cleaned_target} --chr {params.chromosome} --recode --recode-INFO-all --stdout | bgzip -c > {params.path}{params.chromosome}{params.ext} && tabix -p vcf {params.path}{params.chromosome}{params.ext}
        '''

rule phasing_target_RNA_SNP:
    input:
        cleaned_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_CHR{chromosome}.vcf.gz",
        beagle = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/beagle.05May22.33a.jar"
    output:
        phased_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_PHASED_CHR{chromosome}.vcf.gz",
    threads: 40

    params:
        chr = "{chromosome}",
        path = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_PHASED_CHR",
        ext = ".vcf.gz"
    shell:
        """
        java -Xmx110g -jar {input.beagle} gt={input.cleaned_target_chr} impute=false ap=true gp=true nthreads={threads} \
        out={params.path}{params.chr}
        """

#########################################################
# Imputation with variants called from RNA-seq data
#########################################################

rule imputation:
    input:
        minimac4 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/minimac4",
        hap_files = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2.m3vcf.gz",
        phased_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_RNA_PHASED_CHR{chromosome}.vcf.gz"
    output:
        imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz",
        information = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.info"
    threads: 40

    params:
       chr = "{chromosome}",
       prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW"
    shell:
        '''
        {input.minimac4} --refHaps {input.hap_files} --myChromosome {params.chr} --haps {input.phased_target_chr} --prefix {params.prefix}
        '''

rule merge_info_files_RNA:
    input:
        info = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.info", chromosome = autosomes)
    output:
        all_info = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_RNA_ALL_imputation_summary_statistics.info",
        all_sorted = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_RNA_ALL_imputation_summary_statistics_sorted.info"

    shell:
        '''
        cat {input.info} >> {output.all_info}
        sort -n -k1 {output.all_info} > {output.all_sorted}
        '''

rule plot_imputation_performance:
    input:
        RNA_SNPs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_RNA_ALL_imputation_summary_statistics_sorted.info",
        Array_SNPs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/Imputation_performance.R"
    output:
        pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/plots/Imputation_performance_DNA_V_DNA_RNA.pdf"

    shell:
        """
        Rscript {input.script} {input.RNA_SNPs} {input.Array_SNPs} {output.pdf}
        """

rule merge_imputation:
    input:
        imputed = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes)

    output:
        merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.RAW.dose.vcf.gz",

    shell:
        '''
        bcftools concat -Oz -o {output.merged} {input.imputed} && tabix -p vcf {output.merged}
        '''


rule filter_RNA_imputed:
    input:
        renamed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.RAW.dose.vcf.gz"

    output:
        filtered_imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.FILTERED.dose.vcf.gz"
    shell:
        '''
        bcftools +fill-tags {input.renamed} -Oz | bcftools view -q 0.05:minor -Oz | bcftools view -e 'HWE < 0.000001 || R2 < 0.6 || F_MISSING > 0.05 || N_MISSING > 0.05' -Oz -o {output.filtered_imputed} && tabix -p vcf {output.filtered_imputed}
        '''

rule overlap_imputed_variants:
    input:
        filtered_RNA_imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.FILTERED.dose.vcf.gz",
        filtered_DNA_imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz"
    output:
        isec_output = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/sites.txt"
    params:
        file_location = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/"

    shell:
        """
        bcftools isec -p {params.file_location}isec_output -Oz {input.filtered_RNA_imputed} {input.filtered_DNA_imputed}
        """
rule LD_new_and_overlapped:
    input:
        private_RNA = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/0000.vcf.gz", # private to RNA
        ALL_RNA = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_ALL_CHR.IMPUTED.FILTERED.dose.vcf.gz" # overlap
    output:
        private_SNPs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/RNA_private_post.txt",
        ld_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/private_RNA.ld"
    shell:
        """
        bcftools query -f '%ID\n' {input.private_RNA} > {output.private_SNPs}
        
        plink --vcf {input.ALL_RNA} --cow --keep-allele-order --ld-window 15000 --ld-window-kb 1000 --r2 yes-really --ld-window-r2 0.2 --ld-snp-list {output.private_SNPs} --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/isec_output/private_RNA
        """
    



#########################################################
# cis-heritability estimation
#########################################################

rule edit_covariates_gcta:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/Edit_covariates_gcta.R",
        cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt",
        bim = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.bim", 
        fam = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.fam"


    output:
        discrete = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/discrete_covar.covar",
        continuous = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/continuous_covar.qcovar"

    shell:
        """
        Rscript {input.script} {input.bim} {input.fam} {input.cov} {output.discrete} {output.continuous}
        """



rule generate_add_dom_grm:
    input:
        gcta = "/home/workspace/jogrady/my-bin/gcta64/gcta64",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/cis_heritability.R",
        ALL_cis_eVariants = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt",
        gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz",
        bim = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.bim",
        fam = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.fam",
        discrete = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/discrete_covar.covar",
        continuous = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/continuous_covar.qcovar"
    output:
        pheno_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/add_h2/{gene}.phen",
        snp_file_dom  = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/dom_h2/{gene}.cis.SNPlist"

    params:
        gene = "{gene}",
        temp_direc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/add_h2/",
        temp_direc2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/dom_h2/",
        plink_files = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor",
    shell:
        """
        Rscript {input.script} {input.ALL_cis_eVariants} {input.gtf} {params.plink_files} {input.bim} {input.fam} {params.gene} {params.temp_direc} {params.temp_direc2} \
        {input.gcta}
        """
rule plot_heritability:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/plot_heritability.R"
    output:
        pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/Cis-h2-figure.pdf"
    params:
        path = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/add_h2"
    shell:
        """
        Rscript {input.script} {params.path} {output.pdf}
        """
    
#############################################################

# Additive Dominant Model in eQTLs ##########################
#############################################################

rule AD_independent_eQTL:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/Dominance.R",
        data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt",
        data_independent = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_independent_qtl.txt.gz",
        vcf_all = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz",
        expression = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz",
        covariates_all = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt",
        
        data_CON = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt",
        data_independent_CON = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_independent_qtl.txt.gz",
        IFITM3_CON = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.29.txt.gz",
        RGS10_CON = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.26.txt.gz",
        vcf_con = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL_IMPUTED_UPDATED.vcf.gz",
        expression_con = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL.expr_tmm_inv.bed.gz',
        covariates_con = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL.covariates.txt",
        
        data_INFEC = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt",
        data_independent_INFEC = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_independent_qtl.txt.gz",
        IFITM3_INFEC = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.29.txt.gz",
        RGS10_INFEC = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.26.txt.gz",
        vcf_infec = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED_IMPUTED_UPDATED.vcf.gz",
        expression_infec = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED.expr_tmm_inv.bed.gz',
        covariates_infec = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED.covariates.txt",
    output:
        ALL_AD_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/AD_LRT/ALL_AD_results.txt",
        CON_AD_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/AD_LRT/CON_AD_results.txt",
        INFEC_AD_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/AD_LRT/INFEC_AD_results.txt"
    shell:
        """
        Rscript {input.script} {input.data} {input.data_independent} {input.vcf_all} {input.expression} {input.covariates_all} \
        {input.data_CON} {input.data_independent_CON} {input.IFITM3_CON} {input.RGS10_CON} {input.vcf_con} {input.expression_con} {input.covariates_con} \
        {input.data_INFEC} {input.data_independent_INFEC} {input.IFITM3_INFEC} {input.RGS10_INFEC} {input.vcf_infec} {input.expression_infec} {input.covariates_infec} \
        {output.ALL_AD_results} {output.CON_AD_results} {output.INFEC_AD_results}
        """
        
#############################################################

# Gprofiler overlap ##########################
#############################################################
rule Gprofiler_overlap:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/Gprofiler.R",
        DE_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/ALL_results.txt",
        Gprofiler001 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Gprofiler_results.txt"
    
    output:
        Gprofiler005= "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/plots/Gprofiler_ALL_DE_genes_jitter.pdf",
    
    shell:
        """
        Rscript {input.script} {input.DE_results} {input.Gprofiler001} {output.Gprofiler005}
        """


