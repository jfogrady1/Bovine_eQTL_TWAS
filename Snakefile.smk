include: 'snakepit/DeSeq.smk'
include: 'snakepit/Prepare_data_for_eQTL.smk'
include: 'snakepit/Preprocessing.smk'
include: 'snakepit/sample_mismatch.smk'
include: 'snakepit/afc.smk'
include: 'snakepit/residuals.smk'
include: 'snakepit/TWAS.smk'

autosomes = [str(i) for i in range(1,30)] # bovine autosome
rule all:
    input:
        # Preprocessing
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out', individual=config["samples"]),
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt', matrix = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_TPM.txt', matrix = config["groups"]),
        multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Admixture_plot.pdf",
        multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".eigenvec", ".eigenval"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/PLOT_AlleleFrequencies.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN.vcf.gz",
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED.vcf.gz", chromosome = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2.m3vcf.gz", chromosome = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_CHR{chromosome}.vcf.gz", chromosome = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/phased/SNP_Data_PHASED_CHR{chromosome}.vcf.gz", chromosome = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/R2_ER2.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.IMPUTED.RAW.dose.vcf.gz",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.RAW.dose.vcf.gz",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz",
        directory("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/density_raw_counts.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/IBS_histogram_all_animals.pdf",
        # Sample check
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/{sample}.mbv_output.txt", sample = config["samples"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = config["samples"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/Sample_check_T064_MBV.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Volcano_plot.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Volcano_plot_5gPCs.pdf",
        #eQTL analysis
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed.gz", matrix = config["groups"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_tss_file.gtf",
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor.fam", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl.txt.gz", group = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt', group = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.parquet', group = config["groups"], chr = autosomes),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.txt.gz', group = config["groups"], chr = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_independent_qtl.txt.gz", group = config["groups"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_eGene_upset.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_replication_correlation.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.ALL.txt.gz",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/numbers_stats_cis.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.ALL.unzipped.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/violin_plot_pval_thresholds.pdf",
        # Trans
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{group}_genotypes_tensor.vcf.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_deseq_log_norm_counts.bed.gz", matrix = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/{matrix}_AFC.txt", matrix = config["groups"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/Context_specific_eQTL_CON_REAC.pdf",
        # Residuals
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed.gz', group = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_eqtlcovs_qtltools.txt', group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_residualised_expression.txt", group = config["groups"]),
        #TWAS
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas}_GWAS_final.txt", gwas = config["gwas_sets"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}.cis_qtl_pairs_ALL_chr.txt.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_cis_eQTLs.txt", group = config["groups"], gwas_groups = config["gwas_sets"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_{gwas_groups}_TWAS.bed", group = config["groups"], gwas_groups = config["gwas_sets"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/{group}_MeTWAS_models_10_mediators/", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/{group}_{gwas_groups}_TWAS_results.txt", group = config["groups"], gwas_groups = config["gwas_sets"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_TWAS_circos.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/TWAS_intersections.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/CONTROL_volcano_TWAS.pdf"
