import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes


rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz")




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
            --maf_threshold_interaction 0.1 \
            --mode cis_nominal
        '''