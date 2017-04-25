#!/bin/bash

# Job name:
#SBATCH --job-name=merge_delete
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=5G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on error
cp tblastx_result.* $SCRATCH
chkfile tblastx_all.txt
cd $SCRATCH
cat tblastx_result.* > tblastx_res.txt
sort -u  tblastx_res.txt > tblastx_all.txt
cd $SUBMITDIR

rm fastaset*
rm tblastx_result.*


