configfile:
    "config.json"


QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',sample=SAMPLES))


rule all:
    input:
        QC_FILES,
        TRIMMED_FQ

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
