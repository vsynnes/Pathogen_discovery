#!/usr/bin/env python
import numpy as np
import sys
import urllib
import os
import zipfile

#For iterating 2 iterators using the length of the longest
from itertools import izip_longest

#Read nodefile
def read_nodes(myfile):
	while True:
		nodeline = myfile.readline()
		if nodeline =="":
			break
		nodeline = nodeline.split('\t|\t')
		taxid = nodeline[0]
		parent = nodeline[1]
		yield taxid, parent

#Read file containing scientific name
def read_names(myfile):
	for line in myfile:
		line = line.strip('\n').split('\t')
		if line[6] == 'scientific name':
			yield line[0], line[2]
		
	
def main():
	#Downloading taxonomy-file
	print "Downloading taxonomy archive taxdmp.zip ...."
	filenames =set(['nodes.dmp', 'names.dmp', 'merged.dmp'])
	taxfile = urllib.URLopener()
	taxfile.retrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip", "taxdmp.zip")

	pwd=os.getcwd()
	print "Extracting files ...."
	with zipfile.ZipFile('taxdmp.zip', 'r') as inzipfile:
		for name in inzipfile.namelist():
			if name in filenames:
				inzipfile.extract(name,pwd)
	#Deleting zipfile
	os.remove('taxdmp.zip')
       
	#Initializing dictionaries (databases)
        nodes = {}
        names = {}
	merged = {}
	
	#Open nodes and names files
	nodeF = open('nodes.dmp','r')
	nameF = open('names.dmp','r')
	print "Reading taxonomy files ....."

	#Read records of merged taxids	
	with open('merged.dmp','r') as mergeF:
		for line in mergeF:
			line = line.strip('\n').split('\t')
			merged[line[0]] = line[2]

	#Reading nodefile and namefile in parallel
	for (node_taxid, parent_id), (name_taxid, name) in izip_longest(read_nodes(nodeF), read_names(nameF)):
		nodes[node_taxid] = parent_id
		if name_taxid:
			names[name_taxid] = name

  	print "Saving datastructures ..."

    	np.save('nodes.npy',nodes)
        np.save('names.npy',names)
	np.save('merged.npy',merged)

if __name__ == "__main__": main()
