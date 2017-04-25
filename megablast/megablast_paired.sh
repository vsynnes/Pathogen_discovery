#!/bin/bash

# Job name:
#SBATCH --job-name=Megablast_paired_end_reads
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
cp -R blastdb $SCRATCH 
cp ../software/fastq_process.py $SCRATCH

#Assuming that paired end files is present  and copy preprocessed files to SCRATCH
cp ../preprocess/forward_pe.fasta $SCRATCH
cp ../preprocess/forward_pe.fastq $SCRATCH/frw_pe.fastq
cp ../preprocess/reverse_pe.fasta $SCRATCH
cp ../preprocess/reverse_pe.fastq $SCRATCH/rev_pe.fastq

## Make sure the results are copied back to the submit directory (see Work Directory below):
chkfile forward_pe.fastq forward_pe.fasta reverse_pe.fastq reverse_pe.fasta 
chkfile forward_singleton.fastq forward_singleton.fasta reverse_singleton.fastq reverse_singleton.fasta

#Prepare for filter
module load blast+

frw_cnt=0
rev_cnt=0
frw_cnt_s=0
rev_cnt_s=0
# adds arguments in the fastq_process command string based on blast output and the existing singeton files (if any exist)
# and removes duplicates in blast output
add_args(){

        if [ -s forward_removal.txt ]
        then
                python_args="$python_args --discard_file frw_removal.txt"
		sort -u forward_removal.txt > frw_removal.txt
        	cp frw_removal.txt $SUBMITDIR/frw_paired_mega.txt		
                if [ -s reverse_removal.txt ]
                then
                        python_args="$python_args --discard_file2 rev_removal.txt"
			sort -u reverse_removal.txt > rev_removal.txt
                	cp rev_removal.txt $SUBMITDIR/rev_paired_mega.txt
                fi
        elif [ -s reverse_removal.txt ]
        then
                python_args="$python_args --discard_file2 rev_removal.fasta"
		sort -u reverse_removal.txt > rev_removal.txt
              	cp rev_removal.txt $SUBMITDIR/rev_paired_mega.txt
        fi
	append_singletons
}

#Check if singleton files exists, copy the files to scratch, and add arguments to the command string
append_singletons(){
        if [ -s $SUBMITDIR/forward_singleton.fastq ]
        then
		cp $SUBMITDIR/forward_singleton.fastq $SCRATCH
        	cp $SUBMITDIR/forward_singleton.fastq $SCRATCH
		#count sequences
		frw_cnt_s=`expr $(cat forward_singleton.fastq | wc -l) / 4`
                python_args="$python_args --frw_singleOut forward_singleton.fastq forward_singleton.fasta"
        	if [ -s $SUBMITDIR/reverse_singleton.fastq ]
        	then
			cp $SUBMITDIR/reverse_singleton.fastq $SCRATCH
                	cp $SUBMITDIR/reverse_singleton.fastq $SCRATCH
			rev_cnt_s=`expr $(cat reverse_singleton.fastq | wc -l) / 4`
                	python_args="$python_args --rev_singleOut reverse_singleton.fastq reverse_singleton.fasta"
        	fi
        elif [ -s $SUBMITDIR/reverse_singleton.fastq ]
        then
		cp $SUBMITDIR/reverse_singleton.fastq $SCRATCH
                cp $SUBMITDIR/reverse_singleton.fastq $SCRATCH
		rev_cnt_s=`expr $(cat reverse_singleton.fastq | wc -l) / 4`
                python_args="$python_args --rev_singleOut reverse_singleton.fastq reverse_singleton.fasta"
        fi
}

#Get into work directory
cd $SCRATCH

#count number of sequences of input file
frw_cnt=`expr $(cat frw_pe.fastq | wc -l) / 4`
rev_cnt=`expr $(cat rev_pe.fastq | wc -l) / 4`
	
#run blastn on both paired end files
blastn -task megablast  -db blastdb/salmondb -query forward_pe.fasta -evalue 0.0000001 -out forward_removal.txt -outfmt "6 qseqid"
blastn -task megablast  -db blastdb/salmondb -query reverse_pe.fasta -evalue 0.0000001 -out reverse_removal.txt -outfmt "6 qseqid"

#Create command string
python_args="python fastq_process.py frw_pe.fastq --format blast --fastq2In rev_pe.fastq"
add_args

echo $(eval $python_args)
     
#Count number of sequences after subtractions and compute the reduction and print stats
if [ -s forward_pe.fastq ]
then
	new_frw_cnt=`expr $(cat forward_pe.fastq | wc -l) / 4`
	new_frw_cnt_s=`expr $(cat forward_singleton.fastq | wc -l) / 4`

	new_frw_per=$(echo "scale=2;(($frw_cnt - $new_frw_cnt) * 100) / $frw_cnt" | bc)
	new_frw_per_s=$(echo "scale=2;(($new_frw_cnt_s - $frw_cnt_s) * 100) / $frw_cnt_s" | bc)
	echo "After MegaBLAST processing, the dataset contained $new_frw_cnt forward paired end reads, "
	echo "corresponding to $new_frw_per % reduction of the original files."
	echo "The forward singleton file increased by $new_frw_per_s %, and contains $new_frw_cnt_s sequences." 
fi
if [ -s reverse_pe.fastq ]
then
	new_rev_cnt=`expr $(cat reverse_pe.fastq | wc -l) / 4`
	new_rev_cnt_s=`expr $(cat reverse_singleton.fastq | wc -l) / 4`

	new_rev_per=$(echo "scale=2;(($rev_cnt - $new_rev_cnt) * 100) / $rev_cnt" | bc)
	new_rev_per_s=$(echo "scale=2;(($new_rev_cnt_s - $rev_cnt_s) * 100) / $rev_cnt_s" | bc)
	echo "After MegaBLAST processing, the dataset contained $new_rev_cnt reverse paired end reads, "
	echo "corresponding to $new_rev_per % reduction of the original files." 
	echo "The reverse singleton file increased by $new_rev_per_s %, and contains $new_rev_cnt_s sequences."

fi
