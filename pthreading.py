#!/usr/bin/env python2
# Simple Protein Threader
# Written by Arjun Srinivasan

from ptessellation import T3LC,filter_target
from pdb_parser import open_pdb_file, parse_pdb_file
#from fasta_parser import open_fasta_file, parse_fasta_file

def thread_sequence(seq, pdbstruct, replace=True):
	""" Threads a protein sequence through a pdb structure. 

	Arguments:
	seq -- Sequence to be threaded.
	pdbstruct -- Template structure.
	replace -- Not used yet, but will flag the option to replace all gaps in sequences with residues from the original template structure.
	"""
	# The -1 is because resseq begins numbering at 1, while indexing in Python begins at 0
	return [T3LC[seq[x['resseq'] - 1 - (pdbstruct['DBREF']['seqbegin']-pdbstruct['DBREF']['dbseqbegin'])]] for x in pdbstruct["ATOM"] if seq[x['resseq'] - (pdbstruct['DBREF']['seqbegin']-pdbstruct['DBREF']['dbseqbegin'])] != '-' ]


def trim_gaps(seq):
	""" Trims gaps of the primary sequence.

	Arguments:
	seq -- List of sequences. Assumes first element is the primary sequence.
	"""
	return [[y for i,y in enumerate(x) if seq[0][i] != '-'] for x in seq]

# test code below here
#a = parse_pdb_file(open_pdb_file('/home/asriniva/Downloads/3KXU.pdb'), ['ATOM', 'DBREF', 'REMARK'], [465])
#a['ATOM'] = filter_target(a, name=['CA'])
#b = parse_fasta_file(open_fasta_file('/home/asriniva/Downloads/3KXU.fasta.txt'))
#print(b)
#print(thread_sequence(b[0][1], a))
