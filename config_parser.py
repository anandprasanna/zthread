#!/usr/bin/env python3
# Configuration parser
# Written by Arjun Srinivasan

from pdb_parser import *
from fasta_parser import *
from ptessellation import *
from pthreading import *
from output_csv import *

from os import mkdir
import sys

if __name__ == '__main__':
	# Evaluate config file
	exec(open(sys.argv[1], 'r').read())
	if PDB_download:
		pdbfile = download_pdb_file(PDB_id)
	else:
		pdbfile = open_pdb_file(PDB_filename)
	pdbdata = parse_pdb_file(pdbfile, records_parse, records_remarks_parse)
	atoms = filter_target(pdbdata, name=atom_name, span=span, chain=chain)
	verts = tessellate(atoms).vertices
	sim_pot = simplex_potential([[atoms[x]['res'] for x in y] for y in verts]) 
	res_pot = residue_potential(len(atoms), verts, sim_pot)
	print(sum(res_pot))
	if task_thread: 	
		if FASTA_download:
			ffile = download_fasta_file(FASTA_id, FASTA_database)
		else:
			ffile = open_fasta_file(FASTA_filename)
		seq = [x[1] for x in parse_fasta_file(ffile)]
		seq = trim_gaps(seq)
		for z in seq[1:]:
			residues = thread_sequence(z, pdbdata)
			simd = simplex_potential([[residues[x] for x in y] for y in verts]) 
			resd = residue_potential(atoms, verts, simd)
			deltaq = array(resd)-array(res_pot)
			if output_csv:
				mkdir(data_dir+task_name)
				csvfile = open_csv_file(data_dir+task_name+"/output.csv")
#				write_threading_output(csvfile, 

