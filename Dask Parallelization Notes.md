# Dask DataFrames Parallelization

Dask is an python library used for parallelization of code, including pandas using [Dask dataframes](https://saturncloud.io/docs/examples/python/dask/collections/qs-dask-collections-dask-dataframe/)

#dask-dataframes #dask-distributed #dask-mpi

[Another tutorial on Dask dataframe from berkley](https://berkeley-scf.github.io/tutorial-dask-future/python-dask#41-dataframes-pandas)

Dask has [lazy evaluation](https://docs.dask.org/en/stable/user-interfaces.html#laziness-and-computing)
- This means it ONLY computes results when explicitly asked to calling the `df.compute()` function
- Instead it creates a *'task graph'* of how it *would* do things when asked

To use with slurm in sbatch mode is simplest using OpenMPI. Interactive mode (salloc + srun) is trickier because you are only allocated a single node

To install the libraries:

```bash
conda install dask
conda install dask distributed -c conda-forge
conda install dask-mpi -c conda-forge
conda install dask-jobqueue -c conda-forge

# to use dask mpi, also need openmpi
```

### Selecting resources in SLURM to use with Dask

[Dask on HPC Blog](https://blog.dask.org/2019/08/28/dask-on-summit)

#### Definitions

[reference blog post](https://medium.com/analytics-vidhya/how-to-efficiently-parallelize-dask-dataframe-computation-on-a-single-machine-1f10b5b02177#ad77) has nice explanation on threads, processes and how dask utilizes them

 Node: is a computer
 - Each Node =  multiple CPUs
 - Each CPU = mutiple processors (cores)
 - Each Core = multiple threads (only conceptual, not physical)
 - Processes = independent programs using threads to parallelize

**Python cannot use mulitple threads for a single task, therefore Dask relies on processes instead.** This is because of the *Global Interpretor Lock* ( #GIL) that is only bypassed if your data is numeric. For numeric data in arrays or pandas, can use multithreading (multiple cores, 1 process). For all other data, including pure pythonic code, need to use multiprocessing (1 thread, multiple processes).

n_partitions (Dask) >= n_processes (machine) = **ntasks** in slurm

#### Trials

Initially I requested 20 cpus across 5 nodes and 200 GB in total. - that's 100 processes?  
Each worker gets 10 GB  

This is too many - more success with less workers (<20) because each worker gets more memory

`--nodes=5 --ntasks=1 --cpus-per-task=15 --mem=200G` 
- 12.5 GB mem per worker
- 15 min completion

`--nodes=5 --ntasks=15 --cpus-per-task=15 --mem=200G`
- 4.17 GB mem per worker
- lots of high memory use/ unmanaged memory use warnings
- completes task in 15 min

`--nodes=5 --ntasks=15 --cpus-per-task=1 --mem=200G`
- 9.09 GB mem per worker
- least number of errors
- 13 min completion

Additionally, I have no idea whether my tasks are being spread across 5 nodes or not
`--ntasks-per-node` is typically set by scheduler

[on selecting number of CPUs, threads and tasks in slurm](https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch)

### Example sbatch script: 

```bash

#!/bin/bash
#SBATCH --job-name=dask-py           # Job name
#SBATCH --nodes=5                    # Run all processes on a 10 node
#SBATCH --ntasks=15                  # Number of tasks (MPI workers)
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=200G                   # Job memory request
#SBATCH --partition=standard,short,long,fat
#SBATCH --output=dask-py%j.log       # Standard output and error log

#SBATCH --mail-user=varada.khot@uni-jena.de     #email user
#SBATCH --mail-type=BEGIN                       #email when job begins
#SBATCH --mail-type=END                         #email when job ends

echo "Running Dask-MPI"

module load mpi/openmpi/4.1.1

source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh

conda activate working-env  # Activate conda environment 

mpirun --oversubscribe -np 15 python3 filter_diamond_get_clusters.py

'''
# mpirun -np argument needs 2 more workers than what you want for your script
# np rank 0 is for the scheduler
# np rank 1 is for python client process (distribution)
# remaining are workers for your script
'''

```

## How to calculate n workers, memory per worker and chunk size 

Suppose I have a file size 20 GB. 
Read in partition size = 100 MB  
therefore, n_partitions (chunks) = 200  

Here, each partition is brought into RAM, computed and sent back to disk. If I have 4 cores, and total RAM > 4* 100MB, then it will be parallelized.

Number of workers should definitely be < 200 - because you don't want idle workers.
Typically you want workers to do more than just 1 task so workers could be 20-40 in this case (5-10 tasks each)

memory per worker **has to be** >= 100 MB, let's say 500 MB - 1 GB to not run into OOM errors
20 * 0.5 GB

Total memory allocation from slurm required = 10 - 20 GB

[Avoid large partitions](https://docs.dask.org/en/latest/best-practices.html#avoid-very-large-partitions)
- more memory req'd per worker
[Avoid large graphs (many tasks)](https://docs.dask.org/en/latest/best-practices.html#avoid-very-large-graphs)
- overhead per task is 200us - 1ms

[Choosing a good chunksize](https://blog.dask.org/2021/11/02/choosing-dask-chunk-sizes)

[Worker memory management](https://distributed.dask.org/en/stable/memory.html) 
- creating futures
- deleting persisting dataframes/data

[Difference between df.persist() and df.compute()](https://distributed.dask.org/en/stable/memory.html#difference-with-dask-compute)
- df.compute() calls on the execution of the graph
	- don't use with many computations or large data unless that data fits comfortably in local memory
- df.persist() calls on the computation up to a point in the code
	- distributes the computed dask collection across the memory available across all processes (cluster)

[Dask Best Practices](https://docs.dask.org/en/latest/best-practices.html#avoid-very-large-partitions)

[Dask DataFrames Best Practices](https://docs.dask.org/en/stable/dataframe-best-practices.html)


### Example python script dask distributed, dask MPI, dask dataframes :

```python
#! usr/bin/python3

'''
For these large dataframes:
0.5. I used "sort -u" on the command line to get just the sequence IDs of the sseqid,
sort them and deduplicate, instead of reading in full blast output
1. I used the dask library (which parallelizes pandas) to read in my tsv files dd.read_csv instead of pd.read_csv
2. set_index to the columns I wanted to merge on and used df.join() instead of df.merge()
a. can also use df.loc and df.to_dict() + df.map
3. df.persist() intermediate staging of dataframes instead of df.compute() used to locally retrieve data
4. df.persist() stores dataframes into memory so delete dataframes from memory after their use
5. repartition only to a smaller number, not '1'
6. parquet, parquet, parquet... 'to_parquet' is where the computation is happening - multiple parititons are useful
'''

```

```python
# IMPORT LIBRARIES

import pandas as pd
from pandas import DataFrame
import dask.dataframe as dd
from dask.distributed import Client
#from dask_jobqueue import SLURMCluster
from dask_mpi import initialize

```

Initialize the client to work with a cluster - tbh I don't know yet what this does but is important for dask-mpi

```python
initialize()

client = Client() # Connect this local process to remote workers
```

INITIAL READ INS

2 dataframes: a huge one (17GB) + a list of sequence IDs (370K) we want to select from the huge DF

```python
# INITIAL READ INS

# diamond_df = dd.read_csv("./MGYP_v_pharokka_blastp.tsv", sep="\t", header=None,blocksize = '64MB',

# names=['qseqid', 'sseqid','pid','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])

clusters_df = dd.read_csv("/home/ri65fin/MCP_struct/mgnify_databases/mgy_cluster_seqs.tsv", sep="\t", header=None, blocksize = '100 MB',

names=['representative','cluster_seq']).set_index('representative')
print(type(clusters_df))
#print(clusters_df.head(n=10), len(clusters_df))

clusters_df = client.persist(clusters_df)

print("xxxx 1...clusters_df read")

uniq_mgyp_id = dd.read_csv("./mgyp_CR_uniq.list", sep="\t", header=None, names=['CR']).set_index('CR')

print(type(uniq_mgyp_id))
  
# uniq_mgyp_id = diamond_df[['sseqid']].sort_values(by='sseqid', inplace=True).drop_duplicates()

print(uniq_mgyp_id.head(n=10), len(uniq_mgyp_id))
uniq_mgyp_id = uniq_mgyp_id.repartition(npartitions=1)
uniq_mgyp_id = client.persist(uniq_mgyp_id)

print("xxxx 2...uniq_mgyp_id read")

```

Merging the dataframes

`# merged_df = uniq_mgyp_id.merge(clusters_df, how = "left", left_on='sseqid', right_on='representative')`

merging large dataframes is too expensive. A better method is to set_index on the merging columns and use a join on a left index instead

```python
# MERGING DATAFRAMES USING JOIN

# another way to select indices based on the index of a diff column
# merged_df = clusters_df.loc[uniq_mgyp_id.index] 

merged_df = uniq_mgyp_id.join(clusters_df,how='left')
merged_df = client.persist(merged_df)

print("xxxx 3... joining dataframes done")
```

deleting unused dataframes from memory speeds up the compute a lot by freeing up memory

```python
# CLEAN UP DATAFRAMES FROM MEMORY AND REPARTITION

del uniq_mgyp_id
del clusters_df

print("xxxx 3.5... deleted prior dfs from memory")

#print(merged_df.head(n=10), len(merged_df))

# repartition this smaller dataframe
merged_df = merged_df.repartition(npartitions=10).reset_index()
merged_df = client.persist(merged_df)

print("xxxx 4... merged_df repartition to 10 done")
```

write out this intermediate result to parquet

this is by-far the most time intensive step

```python
merged_df.to_parquet('uniq_cluster_seqs_2.parquet')
```

**CONTINUE OPERATIONS ON A SMALLER DATAFRAME**

Could delete the merged_df now and read in pandas dataframe here since the data is now pretty small but I didn't...
 
`# merged_df = pd.read_parquet('uniq_cluster_seqs.parquet')`

```python
# CONTINUE OPERATIONS ON A SMALLER DATAFRAME

merged_df['cluster_seq'] = merged_df['cluster_seq'].str.split(";")
merged_df = client.persist(merged_df)

print("xxxx 4.5 ... split seqids by ;")

# explode rows - 1 row per sequence in cluster
merged_df = merged_df.explode('cluster_seq').reset_index(drop=True) #single column into multiple rows
merged_df = client.persist(merged_df)

print("xxxx 5... dataframe exploded")

print(merged_df.head(n=10), len(merged_df))


# reparition and write to tsv
merged_df = merged_df.repartition(npartitions=1)
merged_df = client.persist(merged_df)

print("xxxx 6... merged_df repartitioned to 1")

merged_df.to_csv('uniq_cluster_seqs_2.tsv', sep='\t', index=False)

```
