#!/usr/bin/env python3
# Configuration parser
# Written by Arjun Srinivasan

from pdb_parser import *
from fasta_parser import *
from ptessellation import *
from pthreading import *

import sys

if __name__ == '__main__':
	# Evaluate config file
	eval(open(sys.argv[1]).read())
	if PDB_download:
		pdbfile = download_pdb_file(PDB_id)
	else:
		pdbfile = open_pdb_file(PDB_filename)
	if FASTA_download:
		ffile = download_fasta_file(FASTA_id, FASTA_database)
	else:
		ffile = open_fasta_file(FASTA_filename)
	

