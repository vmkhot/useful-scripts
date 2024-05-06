# Dask Delayed Parallelization

[Dask delayed](https://tutorial.dask.org/03_dask.delayed.html) parallelization is for other, non-array-non-dataframe, functions including for-loops.

#dask-delayed #dask-distributed #dask-mpi

[Another tutorial on delayed from berkley](https://berkeley-scf.github.io/tutorial-dask-future/python-dask#31-using-a-future-via-delayed)


Example from the Dask tutorial:

```python
import dask

@dask.delayed
def process_file(filename):
    data = read_a_file(filename)
    data = do_a_transformation(data)
    destination = f"results/{filename}"
    write_out_data(data, destination)
    return destination

results = []
for filename in filenames:
    results.append(process_file(filename))

dask.compute(results)
```

pseudo code for multiprocessing fasta files:

```python

# LIBRARIES
import os
import dask
from dask.distributed import Client
from dask_mpi import initialize

# CONNECT TO DISTRIBUTED CLUSTER
initialize()
client = Client()

# READ IN TSV OF SEQUENCE IDS
seqids_dict = {}
'''
1. read in tsv
2. get 2nd column of seqids --> drop duplicates
3. convert to dictionary (index on seqids), key:value = MGYP:1
'''

# READ IN FILENAMES
filenames =[]
for f in os.listdir('/path/to/files'):
	filenames.append(f)

@dask.delayed
def process_file(filename):
	with open # file passed to function as f_handle and output as out:
		for line in f_handle:
			if line startswith ">":
				# split line on whitespace
				# get first element
				# check if first element is in seqids_dict
				# get line and next line
				# write out to "out.fa"

results = []

if __name__ == '__main__':
	
	for filename in filenames:
		results.append(process_file(filename))


results.compute(scheduler='multiprocessing')

```

Actual code:

```python

#! usr/bin/python3
# LIBRARIES
import os, gzip
from pathlib import Path
from itertools import islice
import pandas as pd 
import dask
from dask.distributed import Client
from dask_mpi import initialize

# CONNECT TO DISTRIBUTED CLUSTER
initialize()
client = Client()

# THIS FUNCTION IS PARALLELIZED ACROSS MULTIPLE WORKERS
@dask.delayed
def process_file(filename):
	print(f"processing...{filename}")
	basename = str(filename).split('/')[6].split('.')[0]
	with gzip.open (filename, 'rt') as fh, open (f'./{basename}.fa', 'a+') as out: # file passed to function as f_handle and output as out
		for line in fh:
			if line.startswith(">"):
				line=line.strip('\n')
				mgyp = line.lstrip(">").split()[0] # split line on whitespace and get first element
				# print(mgyp)
				if mgyp in seqid_dict: # check if first element is in seqids_dict
					# print("sequence found...")
					out.write(line +'\n'+''.join(islice(fh,1))) # get line and next line and write out
	return filename

# READ IN TSV OF SEQUENCE IDS
'''
1. read in tsv
2. get 2nd column of seqids --> drop duplicates
3. convert to dictionary (index on seqids), key:value = MGYP:1
'''
seqids_dict = {}
# seqid_df = pd.read_csv('./fake_tsv_for_test.tsv', usecols=['cluster_seq'], sep=r'\s+') 	# for testing
seqid_df = pd.read_csv('../diamond_all_v_all/uniq_cluster_seqs.tsv/0.part', usecols=['cluster_seq'], sep=r'\s+')
# print(seqid_df)
seqid_df.sort_values(by=['cluster_seq']).drop_duplicates(inplace=True)
print(len(seqid_df))
seqid_df['value'] = 1

seqid_dict = dict(zip(seqid_df.cluster_seq, seqid_df.value))
# print(seqid_dict)


# READ IN FILENAMES
filenames =[]

for f in os.listdir('./'):
	if f.endswith(".fa.gz"):
		real_file=os.readlink(f)
		filenames.append(real_file)
print(filenames)

results = []


if __name__ == '__main__':
# HERE WE INVOKE THE PARALLELIZED FUNCTION AND SEND THE FILES TO BE PROCESSED
	for f in filenames:
		results.append(process_file(f))

# DO THE COMPUTATION (BRING BACK A LIST OF THE PROCESSED FILES)
dask.compute(results)
print(results)
print("done processing...")

```

SBATCH script for submitting python script to cluster

```bash

#!/bin/bash

#SBATCH --job-name=dask-py           # Job name
#SBATCH --nodes=1                    # Run all processes on a 10 node
#SBATCH --ntasks=7                  # Number of tasks (MPI workers)
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=30G                   # Job memory request
#SBATCH --partition=standard,short,long,fat
#SBATCH --output=dask-py%j.log       # Standard output and error log

#SBATCH --mail-user=varada.khot@uni-jena.de     #email user
#SBATCH --mail-type=BEGIN                       #email when job begins
#SBATCH --mail-type=END                         #email when job ends

echo "Running Dask-MPI"

module load mpi/openmpi/4.1.1

source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh

conda activate working-env  # Activate conda environment 

mpirun --oversubscribe -np 7 python3 filter_large_fasta_parallel.sh 

```

Runtime:
search 10 sequences in 3 11GB fasta files and write out to fasta in parallel: 5 min
search 3.5 million sequences in 3 11GB fasta files and write them in parallel: 6 min


