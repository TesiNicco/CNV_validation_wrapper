# CNV_validation_wrapper
Wrapper of CNValidatron for validation of PennCNV calls based on genotyping array data.

# Description
This wrapper contains a small adaptation of the CNValidatron package available at https://github.com/SinomeM/CNValidatron_fl/tree/main.  
Please follow all steps described there for installation.  

# Fixes
Parallelization is implemented differently here, as using `BiocParallel` package, the jobs were stalling. Here, we implement a different parallelization strategy, running `{n_cpu}` batches independently of each other in parallel. The total number of samples is evenly split across the defined `{n_cpu}`, and each process is submitted. At the moment, it is based on an interactive shell job, but it can be easily adapted to support `SLURM` submission.  

We also needed to fix another function, which was expecting columns that were not present in the input we used. Note that this does not affect the performances.  

# Installation
As this package will run CNValidatron (https://github.com/SinomeM/CNValidatron_fl/tree/main), make sure that is installed and working.  The `preparation.r` script will also validate the installation of the required R packages, and, in case some of them is missing, will install.  

In general, you will need R (tested on R `4.5.0`) and the following packages:  
`library(BiocManager)`, `library(data.table)`, `library(BiocParallel)`, `library(torch)`, `library(stringr)`, `library(luz)`, `library(CNValidatron)`, `library(parallel)`, `library(argparse)`

# How to run
Clone the repository:
```bash
git clone https://github.com/TesiNicco/CNV_validation_wrapper.git
cd CNV_validation_wrapper/bin
```

The `preparation.r` script will make sure the required packages and input files align with what CNValidatron expects. Briefly, it will:
- check the installation of the required packages, and install those that are missing  
- creates tabix indexed files of the raw intensity data, and adapt the `sample` file  
- adapt the `snp` file  
- adapt the `cnv` file  

The `parallelize_validation.r` will submit the CNValidatron package. Depending on how many `{n_cpu}` are defined, it will create `{n_cpu}` batches each including the same number of samples. The batches will then run in parallel, each using a single thread. By default, each batch will be run as:
`Rscript validation_CNV.R --batch {batch_id} --snps /path/to/snps.txt --cnvs /path/to/cnvs.txt --samples /output_directory/{batch_id}/samples.txt --model /path/to/trained_model.rds --outdir /output_directory/{batch_id}`  
Please refer to https://github.com/SinomeM/CNValidatron_fl/tree/main and https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.621 for a complete description of these inputs, and to download the pre-trained models.  

# Benchmarking performances
| Samples |  CNVs  | CPUs | Time (s) |
|---------|--------|------|----------|
|   1     |   13   |  1   |    14    |
|   5     |   144  |  1   |    88    |
|  10     |   235  |  1   |    144   |
|  45     |   900  |  1   |   3180   |
--------------------------------------
*Table 1: Example benchmarking results for CNV_validation_wrapper using different sample sizes and CNV counts.*
