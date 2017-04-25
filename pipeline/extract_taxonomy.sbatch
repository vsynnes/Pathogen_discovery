#!/bin/bash

# Job name:
#SBATCH --job-name=extract_taxonomy
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=00:10:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=20G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

cp ../post_assembly/tblastx_all.txt $SCRATCH
cp ../software/{names.npy,nodes.npy,merged.npy,extract_taxonomy.py} $SCRATCH

cd $SCRATCH
module load python2
python extract_taxonomy.py tblastx_all.txt lineages.csv
cp lineages.csv $SUBMITDIR/../post_assembly
