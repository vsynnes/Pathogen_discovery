#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import argparse
from itertools import izip_longest

#This generator function takes 2 fastq file as input, 2 id-sets to be removed and yields the paired ends files forward fastq, 
# forward fasta, reverse fastq, reverse fasta, and four similar files for singletons. All eight are only used for paired ends processing
def extract_seq (this_file, rem_set1, rem_set2):
	read_seq = False
	while True:
		#Parse file
		#Search for start of record 		
		this_id, read_seq, this_file = read_ids(this_file, read_seq)
		#call read_seqs and return forward sequence
		this_seq, length, this_file = read_seqs (this_file)
		#read quality
		qual, this_file  = read_qual(this_file, length)
			
		yield add_records(this_id, this_seq, qual, rem_set1, rem_set2, False)

        	if read_seq is False:
			break
		#reset read_seq
		read_seq = False

#reads sequence ids
def read_ids (seq_file, read_seq):
	id_str = ""
	for line in seq_file:
       		if line[0] == '@':
                	read_seq = True
			id_str = line [1:-1]
                        break
	return id_str, read_seq, seq_file

#reads the actual sequence
def read_seqs (seq_file):
	counter, seq_len, seq = 0, 0, ""
	for line in seq_file:
    		counter += len(line) - 1
            	if (line[0] == '+') or (line == ""):
                 	break
		else:
			seq += line
	seq_len = counter - 1
        return  seq, seq_len, seq_file

#Reads quality for fastq files
def read_qual (seq_file, qual_len): 
	counter, vals = 0, ""
	for line in seq_file:
      		counter += len(line) - 1
		vals += line
           	if (counter >= qual_len) or (line == ""):
                	break
	return vals, seq_file
		           
def add_records (seq_id, seq, quality, rem1, rem2, singleton):	
	
	#check if sequence should be removed
	if (seq_id != "" and seq_id in rem1):		
		seq_id, seq, quality = "", "", ""
	#or if counterpart should be removed, mark read as singleton
	elif (seq_id != "" and seq_id[:-1] in rem2):
                singleton = True
        
	return seq_id, seq, quality, singleton

def remove_set (file):
	seq_list = set()
	for line in file:
        	if args.format == 'fasta':
              		#if fasta-file is supplied, read id from fasta fil                              
               		if line[0] == '>' and line[1:-1] not in seq_list:
                        	seq_list.add(line[1:-1])     
		elif line != "" and line[:-1] not in seq_list:#if blast-output is submitted, add ids
               		seq_list.add(line[:-1])
	return seq_list

def open_files(): 
	se_files = []
	#Open paired end output files
	pe_files = [open(pe_names[fileobj], 'w') for fileobj in range(4)]
	#If sequences should be appended to forward singletons
	for fileobj in range(2):
		if args.frw_singleOut:
			se_files.append(args.frw_singleOut[fileobj])
			se_names[fileobj] = args.frw_singleOut[fileobj].name
		else:
			se_files.append(open(se_names[fileobj],'w'))
	#and reverse singletons
	for fileobj in range(2):
		if args.rev_singleOut:
                	se_files.append(args.rev_singleOut[fileobj])
             		se_names[fileobj+2] = args.rev_singleOut[fileobj].name
             	else:
                 	se_files.append(open(se_names[fileobj+2],'w'))
	
	return pe_files, se_files
	

def write_files(id_str, seq, val, fq_file, fa_file):
	#Write fastafile
	if id_str != "":
		fa_file.write('>' + id_str+ '\n')
		fa_file.write(seq)
		#If not fastq->fasta conversion, write fastq file
		if args.format:
			fq_file.write('@' + id_str+ '\n')
			fq_file.write(seq)
			fq_file.write('+\n')
			fq_file.write(val)
	return fq_file, fa_file

#This function closes files and deletes empty files
def close_files(pe_list, se_list):
	for fileobj in range(4):
       		#closing opened files 
		pe_list[fileobj].close()
                se_list[fileobj].close()
		#Delete empty files
		if os.stat(pe_names[fileobj]).st_size==0:
        		os.remove(pe_names[fileobj])
		if os.stat(se_names[fileobj]).st_size==0:
                        os.remove(se_names[fileobj])

	#input files
	if args.fastq2In:
		args.fastq2In.close()
	args.fastqIn.close()

		
if __name__ == "__main__":
	
	
	#Defining commandline arguments--------------------------------------------
	parser = argparse.ArgumentParser()
	parser.add_argument('--format', choices=['blast','fasta'], help='Format of the file containing the sequences to be removed.')
	parser.add_argument('fastqIn',type=argparse.FileType('r'), help='Fastq input file.')
	parser.add_argument('--fastq2In',type=argparse.FileType('r'), help='Second fastq input file if in paired ends mode.')
	parser.add_argument('--orientation',choices=['frw','rev'], help='When processing a single file, indicate if the reads has forward or reverse orientation.')
	parser.add_argument('--min_seq_len', type=int, default=50, help='Minimum length of sequences.')
	parser.add_argument('--discard_file', type=argparse.FileType('r'), help='input file of sequence ids from blast output or fasta sequences to be removed from dataset')
	parser.add_argument('--discard_file2', type=argparse.FileType('r'), help='Second input file of blast output or fasta sequences to be removed from dataset using paired ends.')
	parser.add_argument('--frw_singleOut',nargs=2, type=argparse.FileType('a'), help='The arguments [--frw_single fastq-file fasta-file] loads singleton files in order to append additional sequences to the files.')
	parser.add_argument('--rev_singleOut', type=argparse.FileType('a'), nargs=2, help='The arguments [--rev_single fastq-file fasta-file] loads singleton files in order to append additional sequences to the files.')
	
	#If --discard_file is next argument or no discard_file is indicated and format is, display usage
	if len(sys.argv) < 5: 
		if '--discard_file' in sys.argv or '--format' in sys.argv:
    			parser.print_help()
    			sys.exit(1)
	if '--fastq2In' not in sys.argv and '--orientation' not in sys.argv:
		parser.print_help()
             	sys.exit(1)
	args=parser.parse_args()

	#------------------------------------------------------------------------------------------------#
	
	#Defining filenames to use when writing to file
	se_names = ['forward_singleton.fastq', 'forward_singleton.fasta', 'reverse_singleton.fastq', 'reverse_singleton.fasta']
	pe_names = ['forward_pe.fastq', 'forward_pe.fasta', 'reverse_pe.fastq', 'reverse_pe.fasta']
	
	#sets for storing ids of corresponding sequences to be removed
	remove = set()
	remove2 = set()	
	
	if args.format:
		#Create id-set of sequences to be removed
		remove = remove_set(args.discard_file)

       	 
	pe_files, se_files = open_files()
	if args.format and args.discard_file2:
             	#create second set of unique sequence-id where paire ends records should be removed 
         	remove2 = remove_set(args.discard_file2)
	#If processing only forward file
	if args.orientation=='frw':
		for frw_id, frw_seq, frw_qual, frw_s in extract_seq(args.fastqIn, remove, remove2):
			se_files[0], se_files[1] = write_files(frw_id, frw_seq, frw_qual, se_files[0], se_files[1])
	#If processing only reverse file
	elif args.orientation=='rev':
		for rev_id, rev_seq, rev_qual, rev_s in extract_seq(args.fastqIn, remove, remove2):
			se_files[2], se_files[3] = write_files(rev_id, rev_seq, rev_qual, se_files[2], se_files[3])
	else:			
		#Iterate through generated data by calling extract_seq in parallel with both files
		for (frw_id, frw_seq, frw_qual, frw_s), (rev_id, rev_seq, rev_qual, rev_s) in izip_longest(extract_seq(args.fastqIn, remove, remove2), extract_seq(args.fastq2In, remove2, remove)):
			#When frw_s or rev_s is true, record should be written to singleton file
			if frw_s:
				se_files[0], se_files[1] = write_files(frw_id, frw_seq, frw_qual, se_files[0], se_files[1])
			else:				
				pe_files[0], pe_files[1] = write_files(frw_id, frw_seq, frw_qual, pe_files[0], pe_files[1]) 
			#And for reverse files
			if rev_s:
				se_files[2], se_files[3] = write_files(rev_id, rev_seq, rev_qual, se_files[2], se_files[3])
			else:
				pe_files[2], pe_files[3] = write_files(rev_id, rev_seq, rev_qual, pe_files[2], pe_files[3])	
		
	close_files(pe_files, se_files)
						
