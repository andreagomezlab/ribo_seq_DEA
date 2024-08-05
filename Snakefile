configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001.fastq.gz")

QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',sample=SAMPLES))
RSEM = expand(config['output_dir']+'/rsem/{sample}.genes.results', sample=SAMPLES)

rule all:
    input:
        QC_FILES,
        TRIMMED_FQ,        
        RSEM


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
        nthread = config['nthread']
    output:
        directory(config['output_dir']+'/rsem_bowtie_index_'+config['ref_title'])
    shell:
        r"""
            mkdir -p {params.tmp}; 
            rsem-prepare-reference -p {params.nthread} --gtf {input.gtf} {input.fasta} {output}
        """

rule rsem_calculate_expression:
    input:
        R1=config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',
        R2=config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',
        index = directory(config['output_dir']+'/rsem_bowtie_index_'+config['ref_title'])
    output:
        config['output_dir']+'/rsem/{sample}.genes.results'
    params:
        nthread = config['nthread']
    shell:
        r"""
            rsem-calculate-expression --bowtie2 --num-threads {params.nthread} --paired-end {input.R1} {input.R2} {input.index} {output}
        """


