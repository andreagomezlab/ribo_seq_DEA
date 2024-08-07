configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001.fastq.gz")

QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',sample=SAMPLES))
RSEM = expand(config['output_dir']+'/rsem/{sample}.genes.results', sample=SAMPLES)
BAMS = expand(config['output_dir']+'/rsem/{sample}.genome.bam', sample=SAMPLES)

rule all:
    input:
        QC_FILES,
        TRIMMED_FQ,        
        RSEM,
        BAMS


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
            trim_galore --paired --phred33 --cores 8 {input.R1} {input.R2} -o {params.out_dir} 
        """

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
        config['output_dir']+'/rsem/{sample}.genes.results'    
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
        prefix = config['output_dir']+'/rsem/{sample}',
        index = config['output_dir']+'/rsem_index/rsem_bowtie_index_'+config['ref_title']
    shell:
        r"""
            rsem-tbam2gbam {params.index} {input.bam} {prefix}.genome.bam            
            samtools sort {output} -o {output}.sorted.bam
            samtools index {output}.sorted.bam
            rm {input.bam}
            rm {prefix}.genome.bam            
        """

run RSeQC:
    input:
        bam=config['output_dir']+'/rsem/{sample}.sorted.bam'
    output:
        config['output_dir']+'/RSeQC/{sample}'   
    params:
        gtf=config['genome_gtf'] 
        path=config['rseqc_path']
    log:
        config['output_dir']+'/RSeQC/{sample}.log'
    shell:
        r"""
            {params.path}/bam_stat.py -i {input.bam} > {output}.bam_stat.txt
            
        """

run generate_matrix_genes:
    input:
        rsem_files = RSEM
    output:
        config['output_dir']+'/rsem/genes.results'
    params:
        outdir = config['output_dir']+'/rsem/'
    shell:
        r"""
            rsem-generate-data-matrix {params.outdir} > {output}
        """