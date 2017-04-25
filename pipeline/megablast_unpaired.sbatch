#!/bin/bash

# Job name:
#SBATCH --job-name=Megablast_unpaired_end_reads
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=2:00:00
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
cp -R ../megablast/blastdb $SCRATCH 
cp ../software/fastq_process.py $SCRATCH

#Prepare for MegaBLAST
module load blast+

frw_cnt=0
rev_cnt=0


# check if forward file was returned from  preprocess stage
if [ -s ../preprocess/forward_singleton.fasta ]
then
	#Copy input files
	cp ../preprocess/forward_singleton.fasta $SCRATCH
	#renaming fastq file
	cp ../preprocess/forward_singleton.fastq $SCRATCH/frw_singleton.fastq
	
	#Get into work directory
	cd $SCRATCH
	#count number of sequences of input file
	frw_cnt=`expr $(cat frw_singleton.fastq | wc -l) / 4`
	
	#run megaBLAST
	blastn -task megablast  -db blastdb/salmondb -query forward_singleton.fasta -evalue 0.0000001 -out forward_removal.txt -outfmt "6 qseqid"
	
	#check if removal file exists (blast-result)
	if [ -s forward_removal.txt ]
	then
		#Making file unique
		sort -u forward_removal.txt > frw_removal.txt
		cp frw_removal.txt $SUBMITDIR/../megablast
		#extract host sequences
		python fastq_process.py frw_singleton.fastq --format blast --discard_file frw_removal.txt --orientation frw
		
		#Calculating percentage
		new_frw_cnt=`expr $(cat forward_singleton.fastq | wc -l) / 4`
        	new_frw_per=$(echo "scale=2;(($frw_cnt - $new_frw_cnt) * 100) / $frw_cnt" | bc)
        	echo "After MegaBLAST processing, the dataset contained $new_frw_cnt forward singleton reads, "
        	echo "corresponding to $new_frw_per % reduction of the original files."
	else
		#copy back 
		cp frw_singleton.fastq $SUBMITDIR/../megablast/forward_singleton.fastq
		echo "Processing resulted in no reduction of forward singletons, and contains $frw_cnt sequences."
	fi

	cp forward_singleton.fastq forward_singleton.fasta $SUBMITDIR/../megablast
fi

cd $SUBMITDIR
#similar process with reverse singleton file 
if [ -s ../preprocess/reverse_singleton.fasta ]
then
        #Copy input files
        cp ../preprocess/reverse_singleton.fasta $SCRATCH
        cp ../preprocess/reverse_singleton.fastq $SCRATCH/rev_singleton.fastq
	
	cd $SCRATCH
        #count number of sequences of input file
        rev_cnt=`expr $(cat rev_singleton.fastq | wc -l) / 4`

        #run blastn
        blastn -task megablast  -db blastdb/salmondb -query reverse_singleton.fasta -evalue 0.0000001 -out reverse_removal.txt -outfmt "6 qseqid"
	if [ -s reverse_removal.txt ]
	then
		sort -u reverse_removal.txt > rev_removal.txt
                cp rev_removal.txt $SUBMITDIR/../megablast

		python fastq_process.py rev_singleton.fastq --format blast --discard_file rev_removal.txt --orientation rev
		new_rev_cnt=`expr $(cat reverse_singleton.fastq | wc -l) / 4`
	        new_rev_per=$(echo "scale=2;(($rev_cnt - $new_rev_cnt) * 100) / $rev_cnt" | bc)
	        echo "After MegaBLAST processing, the dataset contained $new_rev_cnt reverse singleton reads, "
        	echo "corresponding to $new_rev_per % reduction of the original files."
	else
		mv rev_singleton.fastq reverse_singleton.fastq
		echo "Processing resulted in no reduction of reverse singletons, and contains $rev_cnt sequences."
	fi
	cp  reverse_singleton.fastq reverse_singleton.fasta $SUBMITDIR/../megablast
fi


























        









        	








