#!/bin/bash

# Job name:
#SBATCH --job-name=merge_delete
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=00:10:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=5G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on error
cp tblastx_result.* $SCRATCH
cd $SCRATCH
cat tblastx_result.* > tblastx_res.txt
sort -u  tblastx_res.txt > tblastx_all.txt
cp tblastx_all.txt $SUBMITDIR/../post_assembly
cd $SUBMITDIR

rm fastaset*
rm tblastx_result.*


