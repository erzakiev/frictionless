# Frictionless (FRIedman transCripTomic sIgnatures Occurring iN muLtiplE Samples Simultaneously)

A blackbox statistical optimization approach to find gene signatures of predefined size in multiple single-cell RNA-seq samples. Computationally demanding, thus better suited for an HPC cluster setting, but a more efficient, single-machine-friendly version of the algorithm is in active development (see [Issue #2](https://github.com/erzakiev/frictionless/issues/2)). A version of the algorithm that doesn't require pre-setting the signature size by the user, using penalization of the objective function, was under active development for a long time, but was largely abandoned due to numerous fruitless attempts of fixing it and due to lack of a relevant background theory in the literature (see [Issue #1](https://github.com/erzakiev/frictionless/issues/1)).

## Installation
### Pulling a [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container
```
singularity pull -U library://erzakiev/default/frictionless:4nov22
```
and then run the `frictionless_4nov22.sif` as an executable:
```
./frictionless_4nov22.sif allInputMatricesRanked 1 Nrestarts k $(date +%s) >> outputfile
```

Manual compilation of the code, downloaded from this repo (the name of the output file is chosen here to be in sync with the containerized version of the pipeline):
```
g++ -O3 main.cpp -o frictionless_4nov22.sif
```

## Usage

## Data format
* Data is expected with cells in columns and genes in rows
* As common practice, you should normalize the reads to individual cell library size, but there is no need to scale the genes
* One operation that needs to be performed by the user, but which will soon be taken care of by the algorithm itself in the nearest releases, is ranking of cells within each gene for each sample separately, with ties split as a mean rank. We recommend R's rank function for that: `input_submatrix_ranked <- t(apply(input_submatrix, 1, rank))`, then concatenating all the submatrices into a single big matrix: `allInputMatricesRanked <- do.call(cbind, list_of_input_submatrices_ranked)`
* Once you have your ranked, concatenated `allInputMatricesRanked`, 
  * Don't forget to set the rownames of it to your appropriate gene names
  * The column names should reflect the cells affiliation to their sample, i.e. cells from the first sample should all have the same name, e.g. `Sample1` or `MGH101`, while cells from the second sample should all have the same, but different from previous sample name, e.g. `Sample2` or `quack`.
    * whatever the chosen notation, cells belonging to the same sample should have the same column name and this is a means of telling the algorithm the affiliation of each cell to a certain sample

We recommend running the algorithm with multiple signature sizes `k`: 10, 20, 30, 40, 50, 100.

### Running on a cluster (recommended)

`run.sh`:
```
#!/bin/bash
## SLURM parameters
##SBATCH -p general # your partition, commented out for a default setting
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --gres=disk:600 # disk space allocated, which should be roughly equal to the size of your input file + 10%
#SBATCH --mem 1200 # memory pool for all cores in MB, which should be roughly equal to twice the size of your input file
#SBATCH -t 0-01:59 # time (D-HH:MM) # 0-00:00 denotes absence of time limit
#SBATCH -o ./slurm.%A.%a.out # STDOUT
#SBATCH -e ./slurm.%A.%a.err # STDERR

## Here is the code that will be run by SLURM
srun  -n $SLURM_NNODES   --ntasks-per-node=1   mkdir /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID} #
sbcast -C allInputMatricesRanked /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID}/allInputMatricesRanked # putting the input file into a temporary local scratch directory for a super fast shared access by all the nodes

srun ./frictionless_4nov22.sif /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID}/allInputMatricesRanked 1 1 k ${{SLURM_JOB_ID}}

srun -n $SLURM_NNODES  --ntasks-per-node=1   rm -rf  /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID} # removing the temp local storage folder


```
and run it with 
```
sbatch --array=0-9999 ./run.sh
``` 
for 10,000 parallel instances.

### Running on a local machine
```
./frictionless_4nov22.sif allInputMatricesRanked 1 Nrestarts k Rand.seed >> outputfile
```
`Nrestarts` is set by the user to the number of restarts your single instance should perform, each restart produces a single signature
`Rand.seed` is a random seed, can be set to the Linux epoch time which is a string of digits `$(date +%s)`: 
```
./frictionless_4nov22.sif allInputMatricesRanked 1 Nrestarts k $(date +%s) >> outputfile
```
