configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001.fastq.gz")

QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
MULTIQC = config['output_dir']+'/qc/multiqc_report.html'
TRIMMED_MULTIQC = config['output_dir']+'/trimmed/multiqc_report.html'
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',sample=SAMPLES))
RSEM_GENES = expand(config['output_dir']+'/rsem/{sample}.genes.results', sample=SAMPLES)
RSEM_ISOFORMS = expand(config['output_dir']+'/rsem/{sample}.isoforms.results', sample=SAMPLES)
BAMS = expand(config['output_dir']+'/rsem/{sample}.sorted.bam', sample=SAMPLES)
COUNTS_G = config['output_dir']+'/rsem/'+config['project_title']+'.genes.results'
COUNTS_I = config['output_dir']+'/rsem/'+config['project_title']+'.isoforms.results'

rule all:
    input:
        QC_FILES,
        #MULTIQC,
        TRIMMED_FQ, 
        #TRIMMED_MULTIQC,       
        RSEM_GENES,
        RSEM_ISOFORMS,
        BAMS,
        COUNTS_G,
        COUNTS_I


rule perform_fastqc:
    input:
        R1=config['fastq']+'/{sample}_R1_001.fastq.gz',
        R2=config['fastq']+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir = config['output_dir']+'/qc'
    output:
        config['output_dir']+'/qc/{sample}_R1_001_fastqc.html',
        config['output_dir']+'/qc/{sample}_R1_001_fastqc.zip',
        config['output_dir']+'/qc/{sample}_R2_001_fastqc.html',
        config['output_dir']+'/qc/{sample}_R2_001_fastqc.zip'
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input.R1} {input.R2}
        '''


rule perform_trim_galore:
    input:
        R1=config['fastq']+'/{sample}_R1_001.fastq.gz',
        R2=config['fastq']+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir = config['output_dir']+'/trimmed'                
    output:
        config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',                    
        config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz'        
    shell:
        r"""
            trim_galore --paired --phred33 --fastqc --cores 8 {input.R1} {input.R2} -o {params.out_dir} 
        """

rule perform_multiqc:
    input:
        expand(config['output_dir']+'/qc/{sample}_R1_001_fastqc.html', sample=SAMPLES),
        expand(config['output_dir']+'/qc/{sample}_R2_001_fastqc.html', sample=SAMPLES)
    output:
        config['output_dir']+'/qc/multiqc_report.html'
    shell:
        r'''
            multiqc -o {config['output_dir']}/qc {input}
        '''

rule rsem_prepare_reference:
    input:
        fasta=config['genome_fasta'],
        gtf=config['genome_gtf']
    params:
        nthread = config['nthread'],
        outdir = config['output_dir']+'/rsem_index/',
        index_prefix = config['output_dir']+'/rsem_index/rsem_bowtie_index_'+config['ref_title']
    output:
        directory(config['output_dir']+'/rsem_index/')        
    shell:
        r"""
            mkdir -p {params.outdir};             
            rsem-prepare-reference -p {params.nthread} --bowtie2 --gtf {input.gtf} {input.fasta} {params.index_prefix}
            
        """

rule rsem_calculate_expression:
    input:
        R1=config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',
        R2=config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',        
        index=config['output_dir']+'/rsem_index/'
    output:
        config['output_dir']+'/rsem/{sample}.genes.results',    
        config['output_dir']+'/rsem/{sample}.isoforms.results',
        config['output_dir']+'/rsem/{sample}.transcript.bam'
    params:
        nthread = config['nthread'],
        prefix = config['output_dir']+'/rsem/{sample}',
        index = config['output_dir']+'/rsem_index/rsem_bowtie_index_'+config['ref_title']
    shell:
        r"""
            rsem-calculate-expression --bowtie2 --num-threads {params.nthread} --paired-end {input.R1} {input.R2} {params.index} {params.prefix}
        """

rule rsem2bam:
    input:        
        bam=config['output_dir']+'/rsem/{sample}.transcript.bam'        
    output:
        config['output_dir']+'/rsem/{sample}.sorted.bam'
    params:        
        index = config['output_dir']+'/rsem_index/rsem_bowtie_index_'+config['ref_title']
    shell:
        r"""
            rsem-tbam2gbam {params.index} {input.bam} {input.bam}.genome.bam
            samtools sort {input.bam}.genome.bam -o {output}
            samtools index {output}
            rm {input.bam}
            rm {input.bam}.genome.bam
        """

rule RSeQC:
    input:
        bam=config['output_dir']+'/rsem/{sample}.sorted.bam'
    output:
        config['output_dir']+'/RSeQC/{sample}'   
    params:
        gtf=config['genome_gtf'], 
        path=config['rseqc_path']
    log:
        config['output_dir']+'/RSeQC/{sample}.log'
    shell:
        r"""
            mkdir -p {output}
            {params.path}/geneBody_coverage.py -r {params.gtf} -i {input.bam} -o {output}/geneBody_coverage
            {params.path}/junction_annotation.py -r {params.gtf} -i {input.bam} -o {output}/junction_annotation
            {params.path}/read_distribution.py -r {params.gtf} -i {input.bam} -o {output}/read_distribution
            {params.path}/RPKM_saturation.py -r {params.gtf} -i {input.bam} -o {output}/RPKM_saturation
        """

rule generate_matrix_genes:
    input:
        rsem_files = RSEM_GENES
    output:
        config['output_dir']+'/rsem/'+config['project_title']+'.genes.results'
    shell:
        r"""
            rsem-generate-data-matrix {input} > {output}
        """

rule generate_matrix_isoforms:
    input:
        rsem_files = RSEM_ISOFORMS
    output:
        config['output_dir']+'/rsem/'+config['project_title']+'.isoforms.results'
    shell:
        r"""
            rsem-generate-data-matrix {input} > {output}
        """