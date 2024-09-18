import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
#rule all:
 #   input:
  #      '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html',
  #      expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out', individual=config["samples"]),
  #      '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
  #      '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
  #      expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt', matrix = config["groups"]),
   #     expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_TPM.txt', matrix = config["groups"]),
  #      multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
 #       "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
 #       "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Admixture_plot.pdf",
    #    multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".eigenvec", ".eigenval"),
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/PLOT_AlleleFrequencies.pdf",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN.vcf.gz",
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED.vcf.gz", chromosome = autosomes),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2.m3vcf.gz", chromosome = autosomes),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_CHR{chromosome}.vcf.gz", chromosome = autosomes),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/phased/SNP_Data_PHASED_CHR{chromosome}.vcf.gz", chromosome = autosomes),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes),
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/R2_ER2.pdf",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.IMPUTED.RAW.dose.vcf.gz",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.RAW.dose.vcf.gz",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz",
    #    directory("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/"),
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/density_raw_counts.pdf",
    #    "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca.pdf",
    #    histogram = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/IBS_histogram_all_animals.pdf"

rule fastqc:
    input:
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/{config["samples"][wildcards.individual]}_{{N}}.fq.gz', N = (1,2))
    output:
        reads = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/'

rule multiqc:
    input:
        reads = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{individual}_{N}_fastqc.zip', individual = config['samples'], N = (1,2))
    output:
        report='/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/
        """
rule build_reference:
    input:
        fa="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf" 
    output:
        STAR_dir = directory("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 149
        ''' 
rule Alignment:
    input:
        genome = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/",
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/{config["samples"][wildcards.individual]}_{{N}}.fq.gz', N=(1,2)),
    params:
        prefix = lambda wildcards: f'{config["samples"][wildcards.individual]}'
    output:
        aligned = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''


rule featureCounts:
    input:
        bam = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam', sample = config["samples"]),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.sh",
        count_matrix_in = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
        final_script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.py",
        temporary = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts_temp.txt"
    output:
        count_matrix = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts2.txt',
        count_matrix_2 = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts_clean.txt',
        cleaned = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
        temporary_1 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts_copy.txt"
    shell:
        '''
        cat {input.count_matrix_in} > {output.temporary_1} 
        tail -n+2 {output.temporary_1} > {output.count_matrix}
        cut -f 1,7-130 {output.count_matrix} > {output.count_matrix_2} # number of samples, get rid of fields we do not want
        sed -i 's#/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/##g' {output.count_matrix_2}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_2} 
        cat {output.count_matrix_2} > {output.cleaned} 
        '''

rule TPM_normalisation:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/TPM_normalisation.R",
        annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        count_matrix_in = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
    output:
        count_matrix_all = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/ALL_matrix_clean.txt',
        count_matrix_control = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/CONTROL_matrix_clean.txt',
        count_matrix_infected = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/INFECTED_matrix_clean.txt',
        count_TPM_all = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/ALL_matrix_TPM.txt',
        count_TPM_control = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/CONTROL_matrix_TPM.txt',
        count_TPM_infected = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/INFECTED_matrix_TPM.txt',
        
    shell:
        '''
        Rscript {input.script} {input.count_matrix_in} {input.annotation} {output.count_matrix_all} {output.count_matrix_control} {output.count_matrix_infected} {output.count_TPM_all} {output.count_TPM_control} {output.count_TPM_infected}
        '''

##### Now onto the Genotype data


rule pruning:
    input:
        plink_files = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data", ".bed", ".bim", ".fam")
    output:
        filtered = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
        pruned = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".bed", ".bim", ".fam")
    shell:
        '''
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data --maf 0.1 --geno 0.95 --mind 0.95 --hwe 0.000001 --allow-extra-chr --double-id --autosome --make-bed --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --indep-pairwise 1000 5 0.2 --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --extract /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset.prune.in --make-bed --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned
        '''

rule admixture:
    input:
       bed_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.bed",
       admixture_exc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/admixture",
    output:
       admixture_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q"
    shell:
        '''
        {input.admixture_exc} --cv {input.bed_file} 2

        mv /home/workspace/jogrady/eqtl_study/eqtl_nextflow/SNP_Pruned* /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/
        '''
rule IBD:
    input:
        pruned = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".bed", ".bim", ".fam"),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Identity_by_state.R"
    output: 
        IBD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.genome",
        histogram = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/IBS_histogram_all_animals.pdf"
    params:
        plink = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned",
    shell:
        '''
        plink --cow --bfile {params.plink} --keep-allele-order --genome --out {params.plink}

        Rscript {input.script} {output.IBD} {output.histogram}
        '''

rule admixture_plot:
    input:
       admixture_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
       plot = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Admixture_plot.R"
    output:
       pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Admixture_plot.pdf"
    shell:
        '''
        Rscript {input.plot} {input.admixture_file} /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE Admixture_plot
        '''

rule pca_plot:
    input:
       filtered = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
       pruned = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset.prune.in",
       admix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
       breed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Sample_database_2022_master.xlsx",
       script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Genotype_PCA_plot.R",
    output:
       pdf = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".eigenvec", ".eigenval"),
       pve = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca_pve.pdf",
       pca = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca.pdf"
    shell:
        '''
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --double-id --allow-extra-chr --extract {input.pruned} --pca --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned 
        Rscript {input.script} {output.pdf[0]} {output.pdf[1]} {output.pve} {input.breed} {input.admix} {output.pca}
        '''

rule liftover:
    input:
       chain = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/bosTau6ToBosTau9.over.chain",
       axiom_master = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Axiom_GW_Bos_SNP_1.na35.annot.csv",
       script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/LiftOver_remapping.R"

    output: 
       liftover = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map"

    shell:
       """
       Rscript {input.script} {input.chain} {input.axiom_master} {output.liftover}
       """

rule remove_palindromic:
    input:
        axiom_master = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Axiom_GW_Bos_SNP_1.na35.annot.csv",
        liftover = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map",
        schnabel = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/9913_ARS1.2_648875_BOS1_marker_name_180910", ".map", ".REF_ALLELE"),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Remove_spurious_SNPs.R"
    output:
        removed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Palindromic_SNPs_removed.txt",
        retained_IDs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Final_raw_IDs_4_analysis.txt",
        liftedover_coord = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_LOver_Rob_final.map",
        final_UMD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt"

    shell:
        '''
        Rscript {input.script} {input.axiom_master} {input.liftover} {input.schnabel} {output.removed} {output.retained_IDs} {output.liftedover_coord} {output.final_UMD}
        '''  
rule remapping:
    input:
        retained_IDs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Final_raw_IDs_4_analysis.txt",
        liftedover_coord = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_LOver_Rob_final.map",
        removed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Palindromic_SNPs_removed.txt",
        schnabel = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/9913_ARS1.2_648875_BOS1_marker_name_180910", ".REF_ALLELE", ".REF"),
        final_UMD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Flipping.R"
    output:
        remapped_plink = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data", ".bed", ".bim", ".fam"),
        remapped_autosomes = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final", ".bed", ".bim", ".fam"),
        vcf_4_flipping = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP.Data.For.Flipping.vcf",
        flip_decision = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ID_Flip_decision.txt",
        final_vcf_4_imputation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        remapped_alleles = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Axiom_Remapped_Flipped_FINAL_A1.vcf",
        flip_coordinates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/regions_file_subset.txt"

    shell:
        '''
        # Update coordinates in plink files

        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data \
        --keep-allele-order --extract {input.retained_IDs} --update-chr {input.liftedover_coord} 2 1 --update-map {input.liftedover_coord} 3 1 \
        --make-bed --allow-extra-chr --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data
        
        # Restrict to autosomes
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data --keep-allele-order \
        --chr 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 --allow-extra-chr --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final

        # Remove genotypes and samples which are missing
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final --keep-allele-order --mind 0.05 --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_1
        
        
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_1 --keep-allele-order --geno 0.05 --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_2
        
        
        # For avoidance of any doubt, set the major allele as reference allele and recode as VCF
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_2 --a1-allele {input.schnabel[0]} --recode vcf \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP.Data.For.Flipping
        
        # Now perform flipping
        Rscript {input.script} {input.final_UMD} {output.vcf_4_flipping} {input.schnabel[1]} {output.flip_decision} {output.final_vcf_4_imputation} {output.remapped_alleles} {output.flip_coordinates}
        
        
        # Redirect files for final 
        cat {output.remapped_alleles} >> {output.final_vcf_4_imputation}
        '''  

rule flip_check:
    input:
        final_vcf_4_imputation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        imputation_master = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED", ".vcf.gz",".vcf.gz.tbi"),
        European_samples = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/European_samples.txt",
        flip_coordinates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/regions_file_subset.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/PLOT_strand_flipping.R"
    output:
        AF_gg = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/PLOT_AlleleFrequencies.pdf",
        allele_frequencies = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/RAW_AF_Comparsion_Values.txt",
        spurious_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Spurious_Flipped_SNPs.txt",
        spurious_snps_pos = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Spurious_SNPS_to_remove_positions.txt"
    shell:
        '''
        # subset genotyped SNPs and European breeds in WGS_reference
        bcftools view -R {input.flip_coordinates} -S {input.European_samples} {input.imputation_master[0]} -Oz \
        -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/WGS_array_European_subsetted.vcf.gz
        
        # use vcftools to calculate frequency
        vcftools --vcf {input.final_vcf_4_imputation} --freq \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Reference_Allele_freq.frq
        
        # Calculate the allele frequencies
        vcftools --gzvcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/WGS_array_European_subsetted.vcf.gz \
        --freq --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/European_WGS_Allele_freq.frq
        
        # Identify spurious SNPs and plot allele frequencies
        Rscript {input.script} \
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/European_WGS_Allele_freq.frq \
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Reference_Allele_freq.frq \
        {output.allele_frequencies} {output.AF_gg} {output.spurious_snps} {output.spurious_snps_pos}
        
        '''

# Final cleaning
rule prepare_ref_target:
    input:
        target_data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        spurious_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Spurious_Flipped_SNPs.txt",
        master_reference = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED", ".vcf.gz",".vcf.gz.tbi")

    params:
        prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/"
    output:
        cleaned_target = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_Final.vcf.gz",
        cleaned_reference = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN", ".vcf.gz", ".vcf.gz.tbi")

    
    shell:
        '''
        # Remove spurious SNPs
        vcftools --vcf {input.target_data} --exclude-positions {input.spurious_snps} --recode --recode-INFO-all --stdout | bgzip -c > {params.prefix}SNP_data_final.vcf.gz
    
        # Removing duplicate SNP - will not work for phasing if present
        gunzip -c {params.prefix}SNP_data_final.vcf.gz | grep -v "Affx-45008270" > {params.prefix}SNP_Data_Final.vcf
    
        gzip {params.prefix}SNP_Data_Final.vcf
    
        # Removing spurious SNPs in WGS reference set
        vcftools --gzvcf {input.master_reference[0]} --exclude-positions {input.spurious_snps} --recode --recode-INFO-all --stdout | bgzip -c > \
        {params.prefix}CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN.vcf.gz \
        && tabix -p vcf {params.prefix}CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN.vcf.gz
        '''

# Now onto the imputation

# Note, WGS reference alread phased from Roslin
rule split_ref_by_chr:
    input:
       cleaned_reference = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED.CLEAN", ".vcf.gz", ".vcf.gz.tbi")
    output:
       imputation_chr = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED", ".vcf.gz", ".vcf.gz.tbi"),
    params:
       chromosome = "{chromosome}"
    shell:
        '''
        vcftools --gzvcf {input.cleaned_reference[0]} --chr {params.chromosome} --recode --recode-INFO-all --stdout | bgzip -c > {output.imputation_chr[0]} && tabix -p vcf {output.imputation_chr[0]}
        '''

rule m3vcfgen:
    input:
       minimac3 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Minimac3",
       imputation_chr = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED", ".vcf.gz", ".vcf.gz.tbi")
    output:
       hap_files = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2.m3vcf.gz"
    
    threads: 40
    params:
       chromosome = "{chromosome}",
       prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2"
    shell:
        '''
        {input.minimac3} --refHaps {input.imputation_chr} --MyChromosome {params.chromosome} --cpus [40]  --processReference --prefix {params.prefix}
        '''

rule split_target_by_chr:
    input:
        cleaned_target = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_Final.vcf.gz"
    output:
        target_chromosome = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_CHR{chromosome}.vcf.gz",
    params:
        chromosome = "{chromosome}",
    shell:
        '''
        vcftools --gzvcf {input.cleaned_target} --chr {params.chromosome} --recode --recode-INFO-all --stdout | bgzip -c > /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_CHR{params.chromosome}.vcf.gz
        '''

rule phasing_target:
    input:
        cleaned_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_CHR{chromosome}.vcf.gz",
        beagle = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/beagle.05May22.33a.jar"
    output:
        phased_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/phased/SNP_Data_PHASED_CHR{chromosome}.vcf.gz"
    threads: 40

    params:
        chr = "{chromosome}"
    shell:
        """
        java -Xmx110g -jar {input.beagle} gt={input.cleaned_target_chr} impute=false ap=true gp=true nthreads={threads} \
        out=/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/phased/SNP_Data_PHASED_CHR{params.chr}
        """

rule imputation:
    input:
        minimac4 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/minimac4",
        hap_files = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/CattleACDr2_noMW_CHR{chromosome}.IMPUTED2.m3vcf.gz",
        phased_target_chr = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/phased/SNP_Data_PHASED_CHR{chromosome}.vcf.gz"
    output:
        imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz"
    threads: 40

    params:
       chr = "{chromosome}",
       prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW"
    shell:
        '''
        {input.minimac4} --refHaps {input.hap_files} --myChromosome {params.chr} --haps {input.phased_target_chr} --prefix {params.prefix}
        '''

# Imputation performance
rule merge_info_files:
    input:
        info = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW.info", chromosome = autosomes)
    output:
        all_info = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics.info",
        all_sorted = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info"

    shell:
        '''
        cat {input.info} >> {output.all_info}
        sort -n -k1 {output.all_info} > {output.all_sorted}
        '''

rule imputation_performance:
    input: 
        all_sorted = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Imputation_Performance.R"
    output:
        ER2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/ER2_values.pdf",
        R2_typed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/R2_typed_variants.pdf",
        ALL_R2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/ALLR2_values.pdf",
        R2_ALL_ER2 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Plotting/imputation/R2_ER2.pdf"

    shell:
        '''
        Rscript {input.script} {input.all_sorted} {output.ER2} {output.R2_typed} {output.ALL_R2} {output.R2_ALL_ER2}
        '''

rule merge_imputation:
    input:
        imputed = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_CHR{chromosome}.IMPUTED.RAW.dose.vcf.gz", chromosome = autosomes)

    output:
        merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.IMPUTED.RAW.dose.vcf.gz",

    shell:
        '''
        bcftools concat -Oz -o {output.merged} {input.imputed} && tabix -p vcf {output.merged}
        '''

rule edit_sample_names:
    input:
        merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.IMPUTED.RAW.dose.vcf.gz",
        names = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Rename.txt"
    
    output:
        renamed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.RAW.dose.vcf.gz",

    shell:
        '''
        bcftools reheader -s {input.names} {input.merged} -o {output.renamed}
        '''
    
rule filter_vcf:
    input:
        renamed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.RAW.dose.vcf.gz"

    output:
        filtered_imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz"
    shell:
        '''
        bcftools +fill-tags {input.renamed} -Oz | bcftools view -q 0.05:minor -Oz | bcftools view -e 'HWE < 0.000001 || R2 < 0.6 || F_MISSING > 0.05 || N_MISSING > 0.05' -Oz -o {output.filtered_imputed} && tabix -p vcf {output.filtered_imputed}
        '''

rule RNA_quality_check:
    input:
        counts = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/RNA_quality_check.R"
    output:
        density = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/density_raw_counts.pdf"

    shell:
        '''
        Rscript {input.script} {input.counts} {output.density}
        '''