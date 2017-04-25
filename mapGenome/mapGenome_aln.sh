#!/bin/bash

# Job name:
#SBATCH --job-name=mapping_Genome
#
# Project:
#SBATCH --account=ln0001k
#
# Wall clock limit:
#SBATCH --time=9:00:00
#

# Number of cores:
#SBATCH --nodes=1
#SBATCH --partition=hugemem
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

## Copy input files to the work directory:
cp ../salmonRef.fa $SCRATCH 
cp ../Sample-2_1.fastq $SCRATCH
cp ../Sample-2_2.fastq $SCRATCH


## Make sure the results are copied back to the submit directory (see Work Directory below):
chkfile mapped_salmon.sorted.bam
chkfile mapped_salmon.sorted.bam.bai
chkfile salmon_unmapped.fasta
chkfile salmon_unmapped.fastq

## Do some work:
cd $SCRATCH

module load bwa
bwa index -a bwtsw salmonRef.fa
bwa aln -t 24 salmonRef.fa Sample-2_1.fastq > alignment1.sai
bwa aln -t 24 salmonRef.fa Sample-2_2.fastq > alignment2.sai
bwa sampe salmonRef.fa alignment1.sai alignment2.sai Sample-2_1.fastq Sample-2_2.fastq > alignment.sam

module load samtools
samtools view -b alignment.sam | samtools sort -o mapped_salmon.sorted.bam
samtools index mapped_salmon.sorted.bam
samtools view  -F 12 -b mapped_salmon.sorted.bam > alignment.bam
samtools fastq alignment.bam > salmon_unmapped.fastq


