#rule all:
 #   input:
  #      expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/{sample}.mbv_output.txt", sample = config["samples"]),
   #     expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = config["samples"]),
    #     "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/Sample_check_T064_MBV.pdf"
rule index_bam:
    input:
       bam = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam",

    output:
        bai = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam.bai"

    shell:
        '''
        samtools index {input.bam}
        ''' 

rule sample_check:
    input:
        vcf = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        check = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''

rule sample_check_plot:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Sample_Check_plot.R",
        check = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/{sample}.mbv_output.txt", sample = config["samples"])
    output:
        pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/Sample_check_T064_MBV.pdf"

    shell:
        '''
        Rscript {input.script} {output.pdf}
        '''