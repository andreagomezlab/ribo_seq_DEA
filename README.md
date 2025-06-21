## System Requirements

To successfully install and run Linux Ubuntu 24.04, your system should have a 2 GHz dual-core processor or better, at least 4 GB of RAM, and a minimum of 25 GB of free hard drive space. The pipeline was tested in a Linux Ubuntu 24.04 AMD Threadripper 3960X 3.8 GHz 24-Core Processor, 256 GB Memory.

## Instructions to use

### Step 1: Copy ribo_seq_DEA workflow repository

#### 1. Clone repository
```shell
$ git clone https://github.com/andreagomezlab/ribo_seq_DEA.git
$ cd ribo_seq_DEA/
```

#### 2. Installing dependencies

#### Conda
Follow the instructions for installing Conda on Linux or other systems

https://docs.anaconda.com/miniconda/miniconda-install/

#### Create and activate the conda environment
```shell
$ cd env/
$ conda env create --name rnaseq_env --file=rnaseq_env.yml 
$ conda activate rnaseq_env
(rnaseq_env)$ pip install --upgrade snakemake
(rnaseq_env)$ snakemake --version
```
Check if the following dependencies are included in the conda environment:
* rsem==v1.3.1
* fastqc==v0.11.9
* trim==galore 0.6.7
* samtools==1.20
* bowtie2==2.4.4
* multiqc==1.14
* snakemake==7.32.4

#### Dry run the pipeline workflow
```shell
$ snakemake -np
```

### Step 2: Configure and execute workflow

Configure the workflow if necessary by editing the file <code>config.json</code> and <code>metadata.cvs</code>

### Execute the workflow

```shell
$(rnaseq_env) snakemake -c10 --use-conda --latency-wait 900

```
### Downstream analysis

RScripts folder
``` shell
edgeR_DEG_analysis_main.R

```

## Demo

### Expected outputs 

```shell
$(rnaseq_env) ls outdir
/rsem/
sample1.sorted.bam
sample1.genes.results
sample1.isoforms.results
/qc/
sample1_R1_001_fastqc.html
sample1_R2_001_fastqc.html
```