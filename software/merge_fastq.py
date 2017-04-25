#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys

import argparse
#This generator function takes 2 fastq file as input, 2 id-sets to be removed and yields the paired ends files forward fastq, 
# forward fasta, reverse fastq, reverse fasta, and four similar files for singletons. All eight are only used for paired ends processing
def extract_seq (frw_file, rev_file):
	read_frw = False
	read_rev = False
	
	#for storing sequence ids
	frw_id = None
	rev_id = None
	

	while True:
		#find first ids in both files		
		if read_frw is False:
			frw_id, read_frw, frw_file = read_ids (frw_file)
		#read reverse file
		if read_rev is False:
			rev_id, read_rev, rev_file = read_ids (rev_file)
		
		#confirm that the fastq files are ordered
		#assert frw_id[:-2] == rev_id[:-2]				
		
		#end of files
		if read_frw is False and read_rev is False:
			break

		#read sequences, call read_seqs and return sequences, length, and file pointer 
		frw_seq, frw_len, frw_file = read_seqs (frw_file)
		rev_seq, rev_len, rev_file = read_seqs (rev_file)
		
		#read quality for fastq file
		frw_qual, frw_file, read_frw = read_qual(frw_file, frw_len)
		rev_qual, rev_file, read_rev = read_qual(rev_file, rev_len)

		#yield id, sequence and quality, interleaving the two files into one
		yield frw_id, frw_seq, frw_qual
                yield rev_id, rev_seq, rev_qual
	
#reads sequence ids
def read_ids (file):
	seq_id, proceed_read = "",  False
	for line in file:
       		if line[0] == '@':
                	proceed_read = True
			fastq_str = line
			seq_id = line
                        break
	return seq_id, proceed_read, file

#reads the actual sequence
def read_seqs (file):
	counter, seq_len, seq = 0, 0, ""
	for line in file:
    		counter += len(line) - 1
		seq += line
            	if (line[:-1] == '+') or (line == ""):
                        break
	seq_len = counter - 1
        return  seq, seq_len, file

#Reads quality for fastq files
def read_qual (file, seq_len): 
	counter, qual, proceed_read = 0, "", False
	for line in file:
      		counter += len(line) - 1
		qual += line
           	if (counter >= seq_len) or (line == ""):    	
                     	break
	return qual, file, proceed_read

		
if __name__ == "__main__":
	
	
	#Defining commandline arguments--------------------------------------------
	parser = argparse.ArgumentParser(description='Script to interleave two fastq input-files')
	parser.add_argument('frw_fastq',type=argparse.FileType('r'), help='Forward input file to be interleaved.')
	parser.add_argument('rev_fastq',type=argparse.FileType('r'), help='Reverse input file to be interleaved.')
	args=parser.parse_args()

	#------------------------------------------------------------------------------------------------#
	#open the file to write sequences to
	out = open('interleaved_seqs.fastq', 'w')
	
	#processing fastq files	by calling generator function extract_seq
	for seq_id, seq, qual in extract_seq(args.frw_fastq, args.rev_fastq):
		#Write id, sequence and quality to file
		out.write(seq_id)
		out.write(seq)
		out.write(qual)
	
	#closing all files
	args.frw_fastq.close()
	args.frw_fastq.close()
	out.close()
						
