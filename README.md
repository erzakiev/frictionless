# Frictionless (FRIedman transCripTomic sIgnatures Occurring iN muLtiplE Samples Simultaneously)

A blackbox statistical optimization approach to find gene signatures of predefined size in multiple single-cell RNA-seq samples. Computationally demanding, thus better suited for an HPC cluster setting, but a more efficient, single-machine-friendly version of the algorithm is in active development (see [Issue #2](https://github.com/erzakiev/frictionless/issues/2)). A version of the algorithm that doesn't require pre-setting the signature size by the user, using penalization of the objective function, was under active development for a long time, but was largely abandoned due to numerous fruitless attempts of fixing it and due to lack of a relevant background theory in the literature (see [Issue #1](https://github.com/erzakiev/frictionless/issues/1)).

## Installation
### Pulling a [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container
```
singularity pull library://erzakiev/default/frictionless
```
and then run the `frictionless_latest.sif` as an executable:
```
./frictionless_latest.sif allInputMatricesRanked Nrestarts k $(date +%s) >> outputfile
```

Manual compilation of the code, downloaded from this repo (the name of the output file is chosen here to be in sync with the containerized version of the pipeline):
```
g++ -O3 main.cpp -o frictionless_latest.sif
```

## Usage

## Data format
* Data is expected with cells in columns and genes in rows
* As per common practice, you should normalize the reads to individual cell library size, but there is *no need* to scale the genes
* One operation that needs to be performed by the user, but which will soon be taken care of by the algorithm itself in the nearest releases, is ranking of cells within each gene for each sample separately, with ties split as a mean rank. We recommend R's rank function for that: `input_submatrix_ranked <- t(apply(input_submatrix, 1, rank))`, then concatenating all the submatrices into a single big matrix: `allInputMatricesRanked <- do.call(cbind, list_of_input_submatrices_ranked)`
* Once you have your ranked, concatenated `allInputMatricesRanked`, 
  * Don't forget to set the rownames of it to your appropriate gene names
  * The column names should reflect the cells affiliation to their sample, i.e. cells from the first sample should all have the same name, e.g. `Sample1` or `MGH101`, while cells from the second sample should all have the same, but different from previous sample name, e.g. `Sample2` or `quack`.
    * whatever the chosen notation, cells belonging to the same sample should have the same column names and this is a means of telling the algorithm which sample each cell came from
    
Here is a dummy usable input matrix `allInputMatricesRanked` with 3 samples called MEF, Repro, Transfo with 10 genes by 10 cells each:
```
MEF MEF MEF MEF MEF MEF MEF MEF MEF MEF Repro Repro Repro Repro Repro Repro Repro Repro Repro Repro Transfo Transfo Transfo Transfo Transfo Transfo Transfo Transfo Transfo Transfo
Dnm3 8 2 5 10 7 4 6 2 2 9 3 9 7 1.5 1.5 10 6 5 8 4 4 2 2 2 6 9 7 8 5 10
Ubald2 7 3.5 9 3.5 3.5 8 3.5 3.5 10 3.5 2 5 6 10 2 2 4 8 7 9 9 2 2 4 10 7 5 8 6 2
Col6a1 4 3 1 8 9 5 2 6 7 10 7 6 3 3 3 3 3 9 8 10 6 2 4 5 10 8 3 1 9 7
Clec4a2 3 8 3 3 6 10 7 3 3 9 3 1.5 10 5 6 4 7 9 1.5 8 4 2 2 8 7 2 10 9 6 5
Gm42372 10 1.5 7 5 4 9 1.5 6 8 3 7 1.5 3 1.5 10 4 8 5 9 6 6 2 5 7 1 3 4 8 10 9
Ret 2 4 10 2 8 9 7 5 2 6 9 5 4 7 10 2 6 3 8 1 8 10 6 7 1 4 3 5 9 2
Sema3d 10 1 7 3 6 8 5 2 4 9 9 2.5 7 2.5 5 6 2.5 8 10 2.5 4 6 9 8 10 7 5 2 2 2
Rian 4 3 2 9 10 6 5 7 8 1 10 2 5 1 6 4 8 7 9 3 9 7 2.5 2.5 6 8 5 2.5 10 2.5
Dsc3 5 9 1 3 10 2 6 7 8 4 10 7 9 2.5 2.5 8 5 2.5 6 2.5 6 3 3 9 10 3 7 3 8 3
Nudt18 6 10 2.5 8 2.5 7 2.5 5 9 2.5 10 1.5 6 8 7 4 9 5 1.5 3 2 10 2 8 4 5 7 9 6 2
```

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

srun ./frictionless_latest.sif /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID}/allInputMatricesRanked Nstarts k ${{SLURM_JOB_ID}}

srun -n $SLURM_NNODES  --ntasks-per-node=1   rm -rf  /local/scratch/tmp/${LOGNAME}_${SLURM_JOB_ID} # removing the temp local storage folder


```
and run it with 
```
sbatch --array=0-9999 ./run.sh
``` 
for 10,000 parallel instances.

### Running one instance on a local Linux machine (ineffecient, not recommended)
```
./frictionless_latest.sif allInputMatricesRanked Nrestarts k $RANDOM >> outputfile
```
`Nrestarts` is set by the user to the number of restarts your single instance should perform, each restart produces a single signature
`$RANDOM` is a random number generated by the OS, seeding the algorithm

### Running on a local Linux machine, using [GNU parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) for parallelization
```
parallel  "./frictionless_latest.sif allInputMatricesRanked Nstarts k {} >> Results_parallel" ::: < RNGs_uniq
```
where we seed each job instance with a pre-prepared file called `RNGs_uniq`, containing random values generated by the `$RANDOM` command as to the total number of jobs desired. For example of a 200-ish jobs (some RNGs repeat, so the final number of unique entries might be slightly smaller than 200) use:
```
for i in {1..200}; do echo $RANDOM >> RNGs; done; cat RNGs | sort --n | uniq > RNGs_uniq
```
A complete, self-contained working example would look like this: 
generating locally using `parallel` ~100 signatures of size 3 genes using 20-ish jobs of 5 restarts each:
```
for i in {1..20}; do echo $RANDOM >> RNGs; done; cat RNGs | sort --n | uniq > RNGs_uniq
parallel  "./frictionless_latest.sif allInputMatricesRanked 5 3 {} >> Results_parallel" ::: < RNGs_uniq
```

And so the output would look like this: each line corresponds to a single signature
```
[signature size] | [total objective function sum] | [indiv. obj. function values for each sample] | [names of genes]
```
```
3	|	4.37	|	1.58	1.15	1.65	|	Col6a1	Rian	Dsc3	
3	|	4.37	|	1.58	1.15	1.65	|	Col6a1	Rian	Dsc3	
3	|	3.91	|	0.77	1.92	1.22	|	Gm42372	Ret	Rian	
3	|	4.37	|	1.58	1.15	1.65	|	Col6a1	Rian	Dsc3	
3	|	3.85	|	1.43	1.37	1.05	|	Rian	Dsc3	Nudt18	
3	|	4.37	|	1.58	1.15	1.65	|	Col6a1	Rian	Dsc3	
3	|	4.37	|	1.58	1.15	1.65	|	Col6a1	Rian	Dsc3	
Program ended with exit code: 0
```


Note: We recommend running the algorithm for the real applications with multiple signature sizes `k`: 10, 20, 30, 40, 50, 100.
