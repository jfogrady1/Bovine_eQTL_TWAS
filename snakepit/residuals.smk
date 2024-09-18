autosomes = [str(i) for i in range(1, 30)]
#rule all:
 #    input:
  #      expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed.gz', group = config["groups"]),
   #     expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_eqtlcovs_qtltools.txt', group = config["groups"]),
    #    expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_residualised_expression.txt", group = config["groups"])



rule edit_bed_file:
    input:
        expression_bed=multiext('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv.bed', '.gz', '.gz.tbi'), # Phenotypes
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/bedfile_edit_qtltools.R" 
    output:
        bed = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed.gz',
        index = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed.gz.tbi'
    params: 
        unzipped = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed',
        group = "{group}"
    shell:
        '''
        Rscript {input.script} {input.expression_bed}[1] {params.group} {params.unzipped} 
        bgzip {params.unzipped}
        tabix -p bed {output.bed}
        '''


rule edit_covariates:
    input:
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Phenotypes
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Covariate_edit_qtltools.R" 
    output:
        cov_edit = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_eqtlcovs_qtltools.txt",
    
    shell:
        '''
        Rscript {input.script} {input.covariates_file} {output.cov_edit} 
        '''


rule extract_residuals:
    input:
        bed = multiext('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/Normalised_counts_{group}_qtltools.bed', ".gz", ".gz.tbi"), # Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_eqtlcovs_qtltools.txt"
    output:
        corrected = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_residualised_expression.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        
        QTLtools correct --bed {input.bed} --cov {input.covariates_file} --out {output.corrected}
        '''
