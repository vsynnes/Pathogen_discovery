#!/usr/bin/env python
import numpy as np
import sys
import os
import time

def main():
#Load taxonomy-databases
    names = np.load('names.npy').item()
    nodes = np.load('nodes.npy').item()
    merged = np.load('merged.npy').item()
    
    infile = sys.argv[1]
    outfile= sys.argv[2]
    taxString = ''
#Open in and out-file
    read = open(infile,'r')
    write = open(outfile,'w')

#Read a line of file and split on 't'
    for line in read:
        line = line.strip('\n').split('\t')
	#In case of multiple taxids, split on ";"
	taxIds = line[1].split(';')
	#Iterate through taxIds
	for thisId in taxIds:
		#Add accession number and taxID to string
        	accession2tax = line[0] + ',' + thisId 
		nextId = thisId
		#if taxid not found
		if nextId not in nodes:
			#Check if the id is merged into another one
			if nextId in merged:
				nextId = merged[nextId]
			else:
				#If taxid does not exist, set taxID to 1 (no rank)
				nextId = '1'
				taxString='no lineage information available'
		#Extract lineage as long as taxID is not 1
		while nextId !='1':
                	taxString = ',' + names[nextId] + taxString
                        nextId = nodes[nextId]
           	#Add lineage behind readId and taxID
        	accession2tax += taxString;
		#write to file
        	write.write(accession2tax+'\n')
       		taxString = ""
    read.close()
    write.close()
if __name__ == "__main__": main()
