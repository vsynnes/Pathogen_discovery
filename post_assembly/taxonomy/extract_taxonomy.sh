#!/bin/bash

# Job name:
#SBATCH --job-name=extract_taxonomy
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=20G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

cp tblastx_all.txt $SCRATCH
cp ../software/{names.npy,nodes.npy,merged.npy,extract_taxonomy.py} $SCRATCH

chkfile lineages.csv

cd $SCRATCH
module load python2
python extract_taxonomy.py tblastx_all.txt lineages.csv

