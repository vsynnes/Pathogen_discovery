#!/bin/bash

# Job name:
#SBATCH --job-name=mapping_Genome
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=14:00:00
#

# Number of cores:
#SBATCH --nodes=1
#SBATCH --partition=hugemem
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

## Copy input files to the work directory:
cp ../salmonref/salmonRef.fa $SCRATCH 
cp ../Sample-2_1.fastq $SCRATCH
cp ../Sample-2_2.fastq $SCRATCH


## Make sure the results are copied back to the submit directory (see Work Directory below):
chkfile left.bam right.bam mapped_salmon.sorted.bam.bai salmon_unmapped1.fastq salmon_unmapped2.fastq alignment.sam

## Do some work:
cd $SCRATCH
raw_cnt_1=`expr $(cat Sample-2_1.fastq | wc -l) / 4`
raw_cnt_2=`expr $(cat Sample-2_2.fastq | wc -l) / 4`
echo "The dataset had $raw_cnt_1 left reads and $raw_cnt_2 right reads prior to processing"

module load bwa
bwa index -a bwtsw salmonRef.fa 
bwa mem  -t 24 salmonRef.fa Sample-2_1.fastq Sample-2_2.fastq > alignment.sam

module load samtools
samtools view -b alignment.sam | samtools sort -o mapped_salmon.sorted.bam
samtools index mapped_salmon.sorted.bam
samtools view  -f 77 -b mapped_salmon.sorted.bam  > left.bam
samtools view  -f 141 -b mapped_salmon.sorted.bam > right.bam
samtools fastq left.bam > salmon_unmapped1.fastq
samtools fastq right.bam > salmon_unmapped2.fastq
left_cnt=`expr $(cat salmon_unmapped1.fastq | wc -l) / 4`
right_cnt=`expr $(cat salmon_unmapped2.fastq | wc -l) / 4`
left_per=$(echo "scale=2;($left_cnt * 100) / $raw_cnt_1" | bc)
right_per=$(echo "scale=2;($right_cnt * 100) / $raw_cnt_2" | bc)
echo "After mapping with BWA the dataset contains $left_cnt left reads and $right_cnt,"
echo " and corresponds to  $left_per % and $right_per % of the complete dataset."
