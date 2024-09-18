## Step 1: Install ribo_seq_DEA workflow

### 1. Clone repository
```shell
$ git clone https://github.com/andreagomezlab/ribo_seq_DEA.git
$ cd ribo_seq_DEA/
```

### 2. Installing Conda

Follow the instructions for installing it on Linux or other systems

https://docs.anaconda.com/miniconda/miniconda-install/


### 3. Create and activate the conda environment
```shell
$ cd env/
$ conda env create --name envname --file=rnaseq_env.yml 
$ conda activate rnaseq_env
```

### 4. Dry run the pipeline workflow
```shell
$ snakemake -np
```

## Step 2: Configure and execute workflow

Configure the workflow if necessary by editing the file <code>config.json</code> and <code>metadata.cvs</code>

## Execute the workflow

```shell
$ snakemake -c10 --use-conda
$ RScripts

```