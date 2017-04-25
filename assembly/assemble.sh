#!/bin/bash

# Job name:
#SBATCH --job-name=BLASTN_paired_end_reads
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=10G

# Number of cores:
#SBATCH --cpus-per-task=1
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

## Copy input files to the work directory
cp ../blastn/forward_pe.fastq $SCRATCH
cp ../blastn/reverse_pe.fastq $SCRATCH
cp ../blastn/forward_singleton.fastq $SCRATCH
cp ../blastn/reverse_singleton.fastq $SCRATCH
cp ../software/merge_fastq.py $SCRATCH

#enter work directory
cd $SCRATCH

#First use flash to merge overlapping sequence pairs
module load flash
#Using the lowest possible maximum overlap value since higher max overlap 
#can reduce  accuracy 
flash --max-overlap=178 forward_pe.fastq reverse_pe.fastq 

#merging singleton files 
python merge_fastq.py forward_singleton.fastq reverse_singleton.fastq
#renaming since running the same sript twice
mv interleaved_seqs.fastq interleaved_singletons.fastq

#merge with extended sequences from flash operation
python merge_fastq.py interleaved_singletons.fastq out.extendedFrags.fastq

#Running velvet by first creating the hashtable using kmer-length from 21 to 49
module load velvet
velveth velvet 23 -fastq -short interleaved_seqs.fastq -fastq -shortPaired -separate out.notCombined_1.fastq out.notCombined_2.fastq
#creating contigs with velvetg
velvetg velvet -min_contig_lgth 50 -exp_cov auto -cov_cutoff 9 -unused_reads yes

#copy back result directory
cp -R velvet $SUBMITDIR
