# System Requirements

# Instructions to use

## Step 1: Copy ribo_seq_DEA workflow repository

### 1. Clone repository
```shell
$ git clone https://github.com/andreagomezlab/ribo_seq_DEA.git
$ cd ribo_seq_DEA/
```

### 2. Installing dependencies

### Conda
Follow the instructions for installing Conda on Linux or other systems

https://docs.anaconda.com/miniconda/miniconda-install/

### Create and activate the conda environment
```shell
$ cd env/
$ conda env create --name rnaseq_env --file=rnaseq_env.yml 
$ conda activate rnaseq_env
(rnaseq_env)$ pip install --upgrade snakemake
(rnaseq_env)$ snakemake --version
```
Check if the following dependencies are included in the conda environment:
* rsem
* fastqc
* trim galore
* samtools
* bowtie2
* multiqc
* snakemake==7.32.4

### Dry run the pipeline workflow
```shell
$ snakemake -np
```

## Step 2: Configure and execute workflow

Configure the workflow if necessary by editing the file <code>config.json</code> and <code>metadata.cvs</code>

## Execute the workflow

```shell
$(rnaseq_env) snakemake -c10 --use-conda --latency-wait 900

```
## Downstream analysis

RScripts folder
``` shell
edgeR_DEG_analysis_main.R

```

## Demo


