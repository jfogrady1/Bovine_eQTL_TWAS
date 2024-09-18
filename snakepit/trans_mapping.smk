import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
#rule all:
    #input:
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal.best.txt.gz", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute.best.txt.gz", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools.bed.gz", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}_trans_FDR_corrected.txt", group = config["groups"]),
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Upset_trans_eGenes.pdf",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/data/links.txt",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot1.png",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_intra.ld.gz",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD.geno.ld",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot2.png",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt"


rule edit_phenotype_bed:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi")

    output:
        bed_4_qtltools = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools", ".bed.gz", ".bed.gz.tbi")

    shell:
        '''
        zcat {input.bed_4_tensor[0]} | awk '{{ $4=$4" . +"; print $0 }}' | tr " " "\t" | bgzip -c > {output.bed_4_qtltools[0]}
        tabix -p bed {output.bed_4_qtltools[0]}
        '''

rule edit_covariates:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Covariate_edit_qtltools.R",
        covs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt"
    output:
        qtl_tools_cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt"
    shell:
        '''
        Rscript {input.script} {input.covs} {output.qtl_tools_cov}
        '''
rule trans_mapping:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz",
        bed_4_qtltools = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools", ".bed.gz", ".bed.gz.tbi"),
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt", # Covariates
    output:
        nominal=multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal",".best.txt.gz", ".bins.txt.gz", ".hits.txt.gz"),
        permute=multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute",".best.txt.gz", ".bins.txt.gz", ".hits.txt.gz")
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       prefix_nominal = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal",
       prefix_permute = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute", 
       group = "{group}"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        # nominal
        QTLtools trans --vcf {input.vcf} --bed {input.bed_4_qtltools[0]} --nominal --cov {input.covariates_file} --normal --threshold 1e-5 --out {params.prefix_nominal}


        # permutation
        QTLtools trans --vcf {input.vcf} --bed {input.bed_4_qtltools[0]} --threshold 1e-5 --cov {input.covariates_file} --normal --permute --out {params.prefix_permute} --seed 1894
        '''
    
rule FDR_correction_trans:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/FDR_trans_correction.R",
        nominal = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal.hits.txt.gz",
        permuted = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute.hits.txt.gz"

    output:
        corrected = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}_trans_FDR_corrected.txt"

    shell:
        '''
        Rscript {input.script} {input.nominal} {input.permuted} {output.corrected}
        '''

rule trans_eGenes_upset:
    input:
        corrected = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}_trans_FDR_corrected.txt", group = config["groups"]),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/trans_upset.R"
    output:
        upset = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Upset_trans_eGenes.pdf"
    shell:
        '''
        Rscript {input.script} {input.corrected} {output.upset}
        '''

rule generate_links:
    input:
        trans_evar = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_trans_FDR_corrected.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Generate_links_broad_trans_eQTL.R"
    output:
        links = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/data/links.txt"

    shell:
        '''
        Rscript {input.script} {input.trans_evar} {output.links}
        '''

rule circos_plot1:
    input:
        config_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/etc/plot1.conf", # note this is made before hand, external to pipeline
        circos = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/bin/circos' # executible
    output:
        plot1_circos = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot1.png"
    shell:
        '''
        {input.circos} -conf {input.config_file}

        mv circos.png /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot1.png
        mv circos.svg /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot1.svg
        '''

rule intrachromosomal_LD_estimate:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Intra_trans_analysis_V2.R"
    output:
        trans_variants = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_chrom_trans_SNPs.txt",
        LD_among_trans_var = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_intra.ld.gz",
        plot_1 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_transVar_LD_between_them.pdf",
        LD_trans_cis_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_observed_cis_trans_SNPs.txt",
        LD_calc_trans_cis_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD.geno.ld",
        null_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/intra_expected_cis_trans_SNPs.txt",
        LD_null = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/trans_cis_intra_LD_null.geno.ld",
        plot_2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_V_15_permuted.pdf",
        permute_Pval = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_cis_permutation_results.txt",
        permuteraw = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/intra_trans_cis_permutation_raw.txt",
        final_intra_trans = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_intrachromosomal_LD_corrected.txt"
    params:
        plink_prefix_path = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor"
    shell:
        '''
        # Note this script performs a few system commands within the script, 
        # in order to calculate LD - please change location of files these if required

        Rscript {input.script} {input.vcf} {input.cis} {input.trans} {output.trans_variants} {output.LD_among_trans_var} \
        {output.plot_1} {output.LD_trans_cis_snps} {output.LD_calc_trans_cis_snps} {output.null_snps} {output.LD_null} \
        {output.plot_2} {output.permute_Pval} {output.permuteraw} {output.final_intra_trans}
        '''
rule interchromosomal_LD_estimate:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/inter_trans_analysis_V2.R"
    output:
        trans_variants_4_LD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs.txt",
        inter_chrom_LD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD.ld",
        null_variants_4_LD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs_NULL.txt",
        null_inter_chrom_LD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD_NULL.ld",
        plot1 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_V_5_permuted.pdf",
        permute_Pval = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_cis_permutation_results.txt",
        permuteraw = "/home/workspace/jograd y/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_cis_permutation_raw.txt",
        final_inter_trans = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt",
        links = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/data/inter_chrom_links.txt"
    shell:
        '''
        # Please check script for actual arguments
        Rscript {input.script}
        '''
rule circos_plot2:
    input:
        config_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/etc/plot2.conf", # note this is made before hand, external to pipeline
        circos = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/bin/circos', # executible,
        links = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/data/inter_chrom_links.txt"
    output:
        plot2_circos = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot2.png"
    shell:
        '''
        {input.circos} -conf {input.config_file}

        mv circos.png /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot2.png
        mv circos.svg /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/circos_plot2.svg
        '''
rule trans_proportion:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/trans_proportion.R",
        inter = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt",
        intra = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_intrachromosomal_LD_corrected.txt",
        TF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_TF.txt",
        coTF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_Cof.txt",
        bed_file ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz",
        vcf_data ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz"
    output:
        trans_prop_stats = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/TF_PROP_STATS.csv",
        trans_prop_raw = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/TF_Prop_RAW.txt",
        pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Trans_TF_proportions.pdf"

    shell:
        """
        # Again, please check inside script
        Rscript {input.script}
        """