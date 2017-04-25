#!/bin/bash

# Job name:
#SBATCH --job-name=blast_contigs
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=10:00:00
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
cp -R ../blastn/blastdb $SCRATCH 
cp ../software/fasta_process.py $SCRATCH

#Assuming that paired end files is present  and copy preprocessed files to SCRATCH
cp ../assembly/velvet/UnusedReads.fa $SCRATCH
cp ../assembly/velvet/contigs.fa $SCRATCH

#Prepare for blast
module load blast+

cd $SCRATCH
#merge files
cat UnusedReads.fa contigs.fa > residual_seqs.fasta

#Count available sequences
cnt=`cat residual_seqs.fasta | grep '>' | wc -l`

blastn -task megablast  -db blastdb/salmondb -query residual_seqs.fasta -evalue 0.0000001 -out mega_removal.txt -outfmt "6 qseqid"

if [ -s mega_removal.txt ]
then
	sort -u mega_removal.txt > mblast_removal.txt
	python fasta_process.py residual_seqs.fasta --format blast --discard_file mblast_removal.txt

	new_cnt=`expr $(cat subtracted_file.fasta | wc -l) / 2`
	new_per=$(echo "scale=2;(($cnt - $new_cnt) * 100) / $cnt" | bc)
	echo "After MegaBLAST processing, the dataset contained $new_cnt reads and contigs, "
	echo "corresponding to $new_per % reduction of the original files."
else
	echo "No sequences was subtracted and the dataset contains $cnt reads and contigs."
fi
#Update old variable
cnt=$new_cnt
#rename outputfile
mv subtracted_file.fasta mega_procs.fasta
blastn -task blastn  -db blastdb/salmondb -query mega_procs.fasta -evalue 0.0000001 -out blast_removal.txt -outfmt "6 qseqid"

if [ -s blast_removal.txt ]
then
	sort -u blast_removal.txt > blastn_removal.txt
        python fasta_process.py mega_procs.fasta --format blast --discard_file blastn_removal.txt

        new_cnt=`expr $(cat subtracted_file.fasta | wc -l) / 2`
        new_per=$(echo "scale=2;(($cnt - $new_cnt) * 100) / $cnt" | bc)
        echo "After BLASTN processing, the dataset contained $new_cnt reads and contigs, "
        echo "corresponding to $new_per % reduction of the original files."
else
        echo "No sequences was subtracted and the dataset contains $cnt reads and contigs."
fi

cp  mblast_removal.txt blastn_removal.txt subtracted_file.fasta $SUBMITDIR/../post_assembly
