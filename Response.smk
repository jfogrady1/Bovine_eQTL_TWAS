import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes


rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz"),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/selected.{sample}.filtered.vcf.gz", sample = config["samples"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/{sample}_rename.txt", sample = config["samples"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_Data_RNA_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/imputation/SNP_RNA_ALL_imputation_summary_statistics_sorted.info"




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
        bcftools reheader -s {output.name_file} {input.out_filter2} -o {output.out_filter3}
        tabix -p vcf {output.out_filter3}
        """

rule merge_RNA_vcf:
    input: 
        vcfs = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/renamed/{sample}.selected.vcf.gz", sample = config["samples"])
    
    output:
        merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_merged.vcf.gz"

    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.merged}
        tabix -p vcf {output.merged}
        """

rule filter_RNA_seq_variant_call:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_merged.vcf.gz"

    output:
        vcf_filtered = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_filtered.vcf.gz"

    shell:
        """
        bcftools view -Oz -m2 -M2 -v snps -i 'F_MISSING < 0.2' -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 {input.vcf} -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
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
        bcftools annotate -x  "INFO/PR" {output.vcf_renamed_header} -Oz -o {output.vcf_renamed_header_noPR}
        tabix -p vcf {output.vcf_renamed_header_noPR}
        """

rule concat_with_SNP_data:
    input:
        RNA_seq = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/RNA_call_filtered.vcf.gz",
        SNP_data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_Data_Final_renamed_noPR.vcf.gz"

    output:
        RNA_SNP_concat = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/merged/SNP_data_RNA_concat_final.vcf.gz"

    shell:
        """
        bcftools concat -a -d snps {input.SNP_data} {input.RNA_seq} -Oz -o {output.RNA_SNP_concat}
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
        vcftools --gzvcf {input.cleaned_target} --chr {params.chromosome} --recode --recode-INFO-all --stdout | bgzip -c > {params.path}{params.chromosome}{params.ext}
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