

#rule all:
 #   input:
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{group}_genotypes_tensor.vcf.gz", group = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_deseq_log_norm_counts.bed.gz", matrix = config["groups"]),
        #expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_AFC.txt", matrix = config["groups"]),
        #context_specific = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/Context_specific_eQTL_CON_REAC.pdf",



rule prepare_genotype_data_afc:
    input:
        vcf = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor.bed")
    output:
        vcf_format = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{group}_genotypes_tensor.vcf.gz"
    
    params:
        plink = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor",
        plink_out = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{group}_genotypes_tensor"
    shell:
        '''
        plink --cow --bfile {params.plink} --keep-allele-order --recode vcf-iid bgz --out {params.plink_out}

        tabix -p vcf {output.vcf_format}      
        '''

rule normalise_deseq2:
    input:
        counts = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt",
        qtl = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{matrix}.cis_qtl.txt.gz", # genes that were only present in qtl mapping
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/afc_normalisation.R",
        covs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.covariates.txt",
        bed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed.gz",

    output:
        norm = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_deseq_log_norm_counts.bed.gz",
        afc_cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}.covariates_qtltools.txt"
    params:
        set = "{matrix}",
        norm_bed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_deseq_log_norm_counts.bed"
    shell:
        '''
        Rscript {input.script} {input.counts} {params.norm_bed} {params.set} {input.qtl} {input.covs} {output.afc_cov} {input.bed}
        bgzip {params.norm_bed}
        tabix -p bed {output.norm}
        '''

rule afc:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_genotypes_tensor.vcf.gz",
        pheno = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_deseq_log_norm_counts.bed.gz",
        cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}.covariates_qtltools.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/aFC.py",
        qtl = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{matrix}.cis_qtl.txt.gz"
    output:
        afc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_AFC.txt"

    shell:
        '''
        python3 {input.script} --vcf {input.vcf} --qtl {input.qtl} --pheno {input.pheno} --cov {input.cov} --log_xform 1 --log_base 2 --geno GT --output {output.afc}
        '''

rule afc_plotting:
    input:
        afc = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_AFC.txt", matrix = config["groups"]),
        eQTL = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{matrix}.cis_qtl_fdr0.05.txt", matrix = config["groups"]),
        nominal_control = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.ALL.unzipped.txt",
        nominal_reactor = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.ALL.unzipped.txt",
        DE = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/ALL_results.txt",
        annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv",
        counts_control = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL_residualised_expression.txt",
        counts_reactor = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED_residualised_expression.txt",
        vcf_control = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL_IMPUTED_UPDATED.vcf.gz",
        vcf_infected = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED_IMPUTED_UPDATED.vcf.gz",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/AFC_eQTL_plotting.r"
    output:
        afc_comparison = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/AFC_comparison_3_groups.pdf",
        afc_maf_slope = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/slope_maf_relationship.pdf",
        afc_slope_cor ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/AFC_slope_correlation.pdf",
        context_specific = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/Context_specific_eQTL_CON_REAC.pdf",
        afc_table_comparison =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/CON_INFEC_comparison.txt",
        IFITM3 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/IFITM3_comparison.pdf",
        IFI16 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/IFI16_comparison.pdf",
        IL1R1 ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/IL1R1_comparison.pdf",
        RGS10 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/RGS10_comparison.pdf"
    shell:
        """
        Rscript {input.script} {input.afc} {input.eQTL} {output.afc_comparison} {output.afc_maf_slope} {output.afc_slope_cor} \
        {input.nominal_control} {input.nominal_reactor} {output.context_specific} {input.DE} {input.annotation} {input.vcf_control} {input.vcf_infected} \
        {input.counts_control} {input.counts_reactor} {output.afc_table_comparison} {output.IFITM3} {output.IFI16} {output.IL1R1} {output.RGS10}
        """