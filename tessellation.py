#!/usr/bin/env python3
# Delaunay Tessellator
# Written by Arjun Srinivasan

from numpy import array
from scipy.spatial import Delaunay
from pdb_parser import parse_pdb_file, download_pdb_file, open_pdb_file

def tessellate(pdbstruct, name='all', span='all', chain='all'):
	""" Performs a Delaunay tessellation on a protein structure.
	
	Arguments:
	pdbstruct -- A dictionary pdb structure as created by the parser.
	name -- A list of atom types to handle. For example, CA.
	span -- Resseq values that are accepted. Values are entered as arguments to range.
	chain -- Chains to parse.
	"""
	# Warnings for missing residues:
	print("The following residues are missing (residue, chain, resseq, icode): "+str(pdbstruct['REMARK'][465]))

	return Delaunay(array([x['coord'] for x in pdbstruct['ATOM'] if (span == 'all' or x['resseq'] in span) and (name == 'all' or x['name'] in name) and (chain == 'all' or x['chain'] in chain)])).vertices


# test code below here
print(tessellate(parse_pdb_file(open_pdb_file('/home/asriniva/Downloads/3KXU.pdb'), ['ATOM', 'REMARK'], [465])))
