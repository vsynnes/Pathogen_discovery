#!/bin/bash

# Job name:
#SBATCH --job-name=trim_reads
#
# Project:
#SBATCH --account =NN9383K
#
# Wall clock limit:
#SBATCH --time=2:00:00
#

# Number of cores:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

## Copy input files to the work directory:
cp ../mapGenome/salmon_unmapped1.fastq $SCRATCH
cp ../mapGenome/salmon_unmapped2.fastq $SCRATCH 
cp ../software/prinseq-lite.pl $SCRATCH
cp ../software/fastq_process.py $SCRATCH 
cp ../software/Trimmomatic-0.36/trimmomatic-0.36.jar $SCRATCH



## Do some work:
cd $SCRATCH

# on which adds arguments in the fastq_process command string
add_args(){
	
        if [ -s bad_paired_1.fasta ]
        then
		#Copy back bad sequences
		cp bad_paired_1.fasta $SUBMITDIR/../preprocess
                python_args="$python_args --discard_file bad_paired_1.fasta"
		if [ -s bad_paired_2.fasta ]
        	then
			#copy back sequences
			cp bad_paired_2.fasta $SUBMITDIR/../preprocess
			python_args="$python_args --discard_file2 bad_paired_2.fasta"
		fi
       	elif [ -s bad_paired_2.fasta ]
	then
		cp bad_paired_2.fasta $SUBMITDIR/../preprocess
		python_args="$python_args --discard_file2 bad_paired_2.fasta"
      	fi
	append_singletons
}

#Helper function which adds singleton arguments in the same manner as above
append_singletons(){
	if [ -s forward_singleton.fastq ]
	then
		python_args="$python_args --frw_singleOut forward_singleton.fastq forward_singleton.fasta"
        	if [ -s reverse_singleton.fastq ]
        	then
        	        python_args="$python_args --rev_singleOut reverse_singleton.fastq reverse_singleton.fasta"
        	fi
	elif [ -s reverse_singleton.fastq ]
       	then
        	python_args="$python_args --rev_singleOut reverse_singleton.fastq reverse_singleton.fasta"
	fi
}

#The code below uses trimmomatic to execute trimming of the paired end reads
#by trimming the ends if lower than phred 3 and uses a sliding window of
#4 nucleotides and cuts read if average phredscore is lower than 5.
java -jar trimmomatic-0.36.jar PE -phred33 salmon_unmapped1.fastq salmon_unmapped2.fastq forward_paired.fastq forward_unpaired.fastq reverse_paired.fastq reverse_unpaired.fastq  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5 MINLEN:50
echo "Checking forward singletons for low complexity reads"

#First forward_unpaired.fastq is produed in  FASTA-format (produces the file forward_singelton.fasta), and then store reads with low
#complexity in bad_forward_unpaired. If this file is not empty, create new 
#updated versions of the original fastq and fasta without the sequences from 
#bad_forward_unpaired. 
if [ -s forward_unpaired.fastq ]
then
	#Convert to fasta
	python fastq_process.py forward_unpaired.fastq --orientation frw
	#Low complexity preprocessing
	./prinseq-lite.pl -fasta forward_singleton.fasta  -out_format 1 -out_bad "bad_forward_unpaired" -lc_method dust -lc_threshold 7

	if [ -s bad_forward_unpaired.fasta ]
	then
		python fastq_process.py --format fasta forward_unpaired.fastq  --discard_file bad_forward_unpaired.fasta --orientation frw
	else
		#first rename forward_unpaired.fastq/a to avoid deletion
		mv forward_unpaired.fastq forward_singleton.fastq
		mv forward_unpaired.fasta forward_singleton.fasta
        	echo "No forward singletons with low complexity preprosessed."
	fi
fi

#Similar type of processing for reverse_unpaired.fastq 
echo "Checking reverse singletons for low complexity reads"

if [ -s reverse_unpaired.fastq ]
then
	#Similar method for reverse_unpaired.fastq, appending to forward singleton file to avoid erasing it
	python fastq_process.py reverse_unpaired.fastq --orientation rev --frw_singleOut forward_singleton.fastq forward_singleton.fasta


	./prinseq-lite.pl -fasta reverse_singleton.fasta  -out_format 1 -out_bad "bad_reverse_unpaired" -lc_method dust -lc_threshold 7
	if [ -s bad_reverse_unpaired.fasta ]
	then
		python fastq_process.py --format fasta reverse_unpaired.fastq  --discard_file bad_reverse_unpaired.fasta --orientation rev --frw_singleOut forward_singleton.fastq forward_singleton.fasta
               
	else
		#first rename revese_unpaired.fastq/a
                mv reverse_unpaired.fastq reverse_singleton.fastq
                mv reverse_unpaired.fasta reverse_singleton.fasta
		echo "No reverse singletons with low complexity preprocessed"
	fi
fi

echo
echo "Checking paired end reads for low complexity reads"
#Processing the set of paired end data under the assumption that 
#data exists in both files

# singleton arguments is necessary to avoid deleting the singleton files made above, since the file name are the same
#and fastq->fasta-conversion yields none singleton sequences 
#In addition, a command string is created and the right type of arguments is inserted for the files present
python_args="python fastq_process.py forward_paired.fastq --fastq2In reverse_paired.fastq" 
append_singletons

#Execute fastq -> fasta conversion
echo $(eval $python_args)

#Extracting low complexity reads
./prinseq-lite.pl -fasta forward_pe.fasta  -fasta2  reverse_pe.fasta -out_format 1 -out_bad "bad_paired" -out_good "paired_proc" -lc_method dust -lc_threshold 7



#Extracting the sequences from the original fastq files
python_args="python fastq_process.py --format fasta forward_paired.fastq --fastq2In reverse_paired.fastq"
add_args
echo $(eval $python_args)
#Copy back files to the preprocess directory
cp forward_pe.fastq forward_pe.fasta reverse_pe.fastq reverse_pe.fasta forward_singleton.fastq forward_singleton.fasta reverse_singleton.fastq reverse_singleton.fasta $SUBMITDIR/../preprocess

