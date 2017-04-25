#!/bin/bash

# Job name:
#SBATCH --job-name=tblastx_nt
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=6:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=10G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

#number of sequences in each file
size="600"
#lines in fasta file
total_lines=`cat ../post_assembly/subtracted_file.fasta | wc -l`
#lines at each process
batch=`echo "$size * 2" | bc`
split -d -l $batch ../post_assembly/subtracted_file.fasta fastaset.
#calculate number of batches, and hence find number of the last process
last_proc=`echo "(($total_lines + $batch - 1) / $batch) - 1" | bc`
arrayrun 00-$last_proc tblastx_remote_worker.sh

