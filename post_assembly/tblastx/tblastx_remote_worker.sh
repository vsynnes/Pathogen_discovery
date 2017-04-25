#!/bin/bash

# Job name:
#SBATCH --job-name=blastn_nt
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=5:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=30G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

FASTASET=fastaset.$TASK_ID
OUTFILE=tblastx_result.$TASK_ID

cp ../pipeline/$FASTASET $SCRATCH
#Copy viruses list
cp viruses.gi $SCRATCH
cd $SCRATCH


#Prepare for blast
module load blast+

cd $SCRATCH

tblastx  -db nt -query $FASTASET -gilist viruses.gi -evalue 0.0001 -max_target_seqs 5 -out $OUTFILE -outfmt "6 qseqid staxids"
cp $OUTFILE $SUBMITDIR/../post_assembly
