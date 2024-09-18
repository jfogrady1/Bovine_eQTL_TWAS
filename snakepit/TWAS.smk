import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
#rule all:
 #   input:
  #      expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas}_GWAS_final.txt", gwas = config["gwas_sets"]),
   #     expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt", group = config["groups"]),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}.cis_qtl_pairs_ALL_chr.txt.gz", group = config["groups"]),
     #   expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_cis_eQTLs.txt", group = config["groups"], gwas_groups = config["gwas_sets"]),
      #  expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_TWAS.bed", group = config["groups"], gwas_groups = config["gwas_sets"]),
       # expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/{group}_MeTWAS_models_10_mediators/", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/{group}_{gwas_groups}_TWAS_results.txt", group = config["groups"], gwas_groups = config["gwas_sets"]),
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_TWAS_circos.pdf",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/TWAS_intersections.txt",
        #"/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/CONTROL_volcano_TWAS.pdf"
rule GWAS_check:
    input:
        ARS_GWAS  = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/{gwas}_ARS.csv",
        UMD_GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/{gwas}_sires.txt",
        Run6 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/1000Bull_UMD_GWAS_REF_ALT_SNPs.tab",
        ARS_alleles = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/1_GWAS_flipping.R"
    output:
        final_GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas}_GWAS_final.txt"
    params:
        set = "{gwas}"      

    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs
        # Note these are massive files, please have the computational resources necessary
        # Otherwise you may need to split by chromosome

        Rscript {input.script} {input.ARS_GWAS} {input.UMD_GWAS} {input.Run6} {input.ARS_alleles} {params.set} {output.final_GWAS}
        '''
rule MeQTL_correlation:
    input:
        exp = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv.bed.gz",
        TF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_TF.txt",
        coTF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_Cof.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/2_MediatorQTLAnalysis.R"

    output:
        med_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_Mediators_location.txt",
        med_intensity = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_Mediators_intensity.txt",
        med_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt"

    shell:
        '''
        Rscript {input.script} {input.exp} {input.TF} {input.coTF} {output.med_loc} {output.med_intensity} {output.med_results}
        '''

rule match_GWAS:
    input:
        thresholds = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_cis_qtl_pairs_TWAS_input.txt",
        correlations = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt",
        GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas_groups}_GWAS_final.txt",
        expression = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv.bed.gz",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/3_match_GWAS_QTL.R"
    output:
        eqtls = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_cis_eQTLs.txt",
        snp_pos = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_snp_pos.txt",
        final_correlations = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_correlations.txt",
        GWAS_out = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_GWAS.txt",
        mediators_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_Mediator_LOCs.txt",
        mediators_exp = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_Mediator_intensities.txt",
        gene_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_genes_LOCs.txt",
        gene_exp = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_genes_intensity.txt",
        snpid = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_final_SNPs.txt"

    shell:
        '''
        Rscript {input.script} {input.thresholds} \
        {input.correlations} \
        {input.GWAS} \
        {input.expression} \
        {output.eqtls} {output.snp_pos} {output.final_correlations} \
        {output.GWAS_out} {output.mediators_loc} {output.mediators_exp} \
        {output.gene_loc} {output.gene_exp} {output.snpid}
        '''

rule makebedfile:
    input:
        snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_final_SNPs.txt",
    output:
        twas_inputs = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_TWAS", ".bed", ".bim", ".fam", ".vcf")
    params:
        plinkfile = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor",
        plink_out = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_TWAS"
    shell:
        '''
        plink --cow --bfile {params.plinkfile} --keep-allele-order --extract {input.snps} --make-bed --out  {params.plink_out}
        plink --cow --bfile {params.plinkfile} --keep-allele-order --extract {input.snps} --recode vcf --out  {params.plink_out}
        '''

rule MeTWAS_generation:
    # Only need to do it once for each group
    input:
        bed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_TWAS.bed",
        med = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_Mediator_intensities.txt",
        exp = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_genes_intensity.txt",
        med_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_Mediator_LOCs.txt",
        gene_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_genes_LOCs.txt",
        covariates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt",
        qtl_local = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_cis_eQTLs.txt",
        correlations = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_ALL_correlations.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/3_MeTWAS_Generation.R"
    output:
        models_all = directory("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/{group}_MeTWAS_models_10_mediators/")
    shell:
        '''
        Rscript {input.script} {input.bed} {input.med} \
        {input.exp} {input.med_loc} {input.gene_loc} {input.covariates} \
        {input.qtl_local} {input.correlations} \
        '''
rule TWAS:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/4_BurdenTest.R",
        GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_GWAS.txt",
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_TWAS.vcf",
        model = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/{group}_MeTWAS_models_10_mediators/",
    output:
        TWAS_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/{group}_{gwas_groups}_TWAS_results.txt"

    shell:
        '''
        Rscript {input.script} {input.GWAS} {input.vcf} {input.model} {output.TWAS_results}
        '''
rule circos:
    input:
        CH = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_CH_TWAS_results.txt",
        LM = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_LM_TWAS_results.txt",
        HF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_HOFR_TWAS_results.txt",
        ALL = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_ALL_TWAS_results.txt",
        ensemble = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/Circos_ALL_TWAS.R"

    output:
        pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_TWAS_circos.pdf"

    shell:
        '''
        Rscript {input.script} {input.CH} {input.LM} {input.HF} {input.ALL} {input.ensemble} {output.pdf}
        '''

rule Upset_TWAS:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/Upset_TWAS.R",
        ALL_CH = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_CH_TWAS_results.txt",
        ALL_LM = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_LM_TWAS_results.txt",
        ALL_HF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_HOFR_TWAS_results.txt",
        ALL_ALL = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_ALL_TWAS_results.txt",
        CON_CH = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_CH_TWAS_results.txt",
        CON_LM = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_LM_TWAS_results.txt",
        CON_HF ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_HOFR_TWAS_results.txt",
        CON_ALL ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_ALL_TWAS_results.txt",
        REA_CH = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_CH_TWAS_results.txt",
        REA_LM = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_LM_TWAS_results.txt",
        REA_HF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_HOFR_TWAS_results.txt",
        REA_ALL = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_ALL_TWAS_results.txt",
        annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv"
    output:
        upset = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/TWAS_upset.pdf",
        intersections = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/TWAS_intersections.txt"
    
    shell:
        '''
        Rscript {input.script} {input.ALL_CH} {input.ALL_LM} {input.ALL_HF} {input.ALL_ALL} \
        {input.CON_CH} {input.CON_LM} {input.CON_HF} {input.CON_ALL} {input.REA_CH} {input.REA_LM} {input.REA_HF} {input.REA_ALL} \
        {output.upset} {input.annotation} {output.intersections}
        '''
rule volcano_TWAS:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/6_TWAS_annotating_plotting.R",
        in_1 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_CH_TWAS_results.txt",
        in_2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_LM_TWAS_results.txt",
        in_3 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_HOFR_TWAS_results.txt",
        in_4 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_ALL_TWAS_results.txt",
        in_5 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_CH_TWAS_results.txt",
        in_6 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_LM_TWAS_results.txt",
        in_7 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_HOFR_TWAS_results.txt",
        in_8 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_ALL_TWAS_results.txt",
        in_9 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_CH_TWAS_results.txt",
        in_10 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_LM_TWAS_results.txt",
        in_11 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_HOFR_TWAS_results.txt",
        in_12 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_ALL_TWAS_results.txt",
        in_13 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv"

    output:
        out_14 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_Supplementary.txt",
        out_15 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_Supplementary.txt",
        out_16 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_Supplementary.txt",
        out_17 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/INFECTED_volcano_TWAS.pdf",
        out_18 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/CONTROL_volcano_TWAS.pdf"

    shell:
        '''
        Rscript {input.script} {input.in_1} {input.in_2} {input.in_3} {input.in_4} {input.in_5} {input.in_6} \
        {input.in_7} {input.in_8} {input.in_9} {input.in_10} {input.in_11} {input.in_12} {input.in_13} \
        {output.out_14} {output.out_15} {output.out_16} {output.out_17} {output.out_18}
        '''