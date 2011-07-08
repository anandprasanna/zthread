#!/usr/bin/env python2
# FASTA file parser
# Written by Arjun Srinivasan

#from urllib.request import urlopen
from urllib import urlopen

UNIPROT_URL = "http://www.uniprot.org/uniprot/"

def open_fasta_file(ffilename):
	""" Opens local FASTA file for reading.

	Arguments:
	ffilename -- FASTA filename.
	"""
	return open(ffilename, 'rb')

def download_fasta_file(fid, database):
	""" Downloads FASTA files from databases (such as UNIPROT or NCBI).

	Arguments:
	fid -- The id of the protein in each respective database.
	database -- The string identifier for the database.
	"""
	if database == 'unp':
		return urlopen(UNIPROT_URL + fid + ".fasta")
	elif database == 'gi':
		print("Not supported because NCBI makes it annoying.")

def parse_fasta_file(ffile):
	""" Parses FASTA files.

	Arguments:
	ffile -- Open filehandle of a fasta file.
	"""
	sequences = str(ffile.read(), 'ascii').split('>')
	ffile.close()
	return [(x.splitlines()[0], ''.join(x.splitlines()[1:])) for x in sequences[1:]]

# test code below here

#print(parse_fasta_file(open_fasta_file("/home/asriniva/Downloads/hemoglobin.fasta.txt")))


