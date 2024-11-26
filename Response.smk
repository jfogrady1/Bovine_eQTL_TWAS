import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes


rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/SplitNCigarReads_C002-STARAligned.sortedByCoord.out.bam"




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
sample = ""  # sample id

# 1. Add read groups, sort, mark duplicates, and create index
rule task1:
    input:
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/C002_Aligned.sortedByCoord.out.bam"
    output:
        read_group = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/rg_added_C002-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/dedupped_C002-STARAligned.sortedByCoord.out.bam"
    threads: 10
    params: TMPDIR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/"
    shell:
        '''        
        java -jar -Xmx15g -Djava.io.tmpdir={params.TMPDIR} {input.picard} AddOrReplaceReadGroups I={input.in_bam} O={output.read_group} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        java -jar -Xmx15g -Djava.io.tmpdir={params.TMPDIR} {input.picard} MarkDuplicates I={output.read_group} O={output.dedupped} \
        CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=marked_dup_metrics.txt
        '''

# 2. Split N trim and reassign mapping qualities
rule task2:
    input:
        reference="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        dedupped ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/dedupped_C002-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/SplitNCigarReads_C002-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/"
    threads: 10
    shell:
        '''
        {input.gatk} SplitNCigarReads -R {input.reference} -I {input.dedupped} -O {output.split} 
        '''
