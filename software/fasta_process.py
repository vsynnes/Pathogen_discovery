#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import argparse

#This generator function takes a fasta file as input, a id-set of sequences to be removed and yields the reduced fasta file. 
def extract_seq (fasta_file, rem_set):
	next_id = ""
	read_rec = False
	while True:
		#Read a record 		
		seq_id, next_id, sequence, fasta_file, read_rec = read_record (fasta_file, next_id)
		#Add record only if it should be retained
		if seq_id not in rem_set:
			yield seq_id, sequence
		if read_rec:
			read_rec = False
		else:
			break

#reads sequence ids
def read_record (seq_file, next_id):
	id_str, seq, read_next = "", "", False
	if next_id == "":
		for line in seq_file:
       			if line[0] == '>':
				id_str = line [1:-1]
                        break
	else:
		id_str = next_id

	for line in seq_file:
		if line[0] == '>':
			next_id = line[1:-1]
			read_next = True
			break
		#If line is not empty string or endline
		elif line[:-1] !="" or line !="":
			seq += line.strip('\n')
		
	return id_str, next_id, seq, seq_file, read_next

		           

def remove_set (file):
	seq_list = set()
	for line in file:
		seq_id = line.strip('\n').split('\t')[0]
        	if args.format == 'fasta':
              		#if fasta-file is supplied, read id from fasta fil                              
               		if seq_id[0] == '>' and seq_id[1:-1] not in seq_list:
                        	seq_list.add(line[1:-1])     
		elif seq_id != "" and seq_id[:-1] not in seq_list:#if blast-output is submitted, add ids
               		seq_list.add(line[:-1])
	return seq_list


		
if __name__ == "__main__":
	
	
	#Defining commandline arguments--------------------------------------------
	parser = argparse.ArgumentParser()
	parser.add_argument('--format', choices=['blast','fasta'], help='Format of the file containing the sequences to be removed.')
	parser.add_argument('fastaIn',type=argparse.FileType('r'), help='Fasta input file.')
	parser.add_argument('--discard_file', type=argparse.FileType('r'), help='input file of sequence ids from blast output or fasta sequences to be removed from dataset')
	
	#If --discard_file is next argument or no discard_file is indicated and format is, display usage
	if len(sys.argv) < 5: 
		if '--discard_file' in sys.argv or '--format' in sys.argv:
    			parser.print_help()
    			sys.exit(1)
	args=parser.parse_args()

	#------------------------------------------------------------------------------------------------#
	
	#sets for storing ids of corresponding sequences to be removed
	remove = set()
	#open output file
	out = open('subtracted_file.fasta','w')
	if args.format:
		#Create id-set of sequences to be removed
		remove = remove_set(args.discard_file)
			
	#Iterate through generated data by calling extract_seq
	for fasta_id, fasta_seq in extract_seq(args.fastaIn, remove):
		#Write sequence string first, then actual sequence
		if fasta_id!="" and fasta_seq!="":
			out.write('>'+fasta_id+'\n')
			out.write(fasta_seq+'\n')
	out.close()
	args.fastaIn.close()				
