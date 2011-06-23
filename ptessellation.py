#!/usr/bin/env python3
# Delaunay Tessellator
# Written by Arjun Srinivasan

from numpy import array
from scipy.spatial import Delaunay
from pdb_parser import parse_pdb_file, download_pdb_file, open_pdb_file

POTENTIAL_DICT_FILE = "potential_1417.out"
TO_ONE_LETTER_CODE={'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M','ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K','ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H','CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'}
TOLC = TO_ONE_LETTER_CODE
print('asdf')
T3LC = {y:x for x,y in TOLC.items()}

def generate_potential_dictionary():
	""" Generate the potential dictionary for four-body potentials based on previous data. """
	with open(POTENTIAL_DICT_FILE, 'r') as f:
		return {x.split()[0] : float(x.split()[1]) for x in f}

POTENTIAL_DICT = generate_potential_dictionary()

def filter_target(pdbstruct, name='all', span='all', chain='all'):
	""" Filters a protein structure based on given information.
	
	Arguments:
	pdbstruct -- A dictionary pdb structure as created by the parser.
	name -- A list of atom types to handle. For example, CA.
	span -- Resseq values that are accepted. Values are entered as arguments to range.
	chain -- Chains to parse.
	"""
	# Warnings for missing residues:
	print("The following residues are missing (residue, chain, resseq, icode): "+str(pdbstruct['REMARK'][465]))

	return [x for x in pdbstruct['ATOM'] if (span == 'all' or x['resseq'] in span) and (name == 'all' or x['name'] in name) and (chain == 'all' or x['chain'] in chain)]


def tessellate(atoms):
	""" Performs a Delaunay tessellation on a given list of atoms.
	
	Arguments:
	atoms -- List of atoms to triangulate.

	Returns:
	A SciPy Delaunay object (contains info about vertices and neighbors)
	"""
	return Delaunay(array([x['coord'] for x in atoms]))


def potential(atomlist):
	return [(''.join([TOLC[y['res']] for y in x]), POTENTIAL_DICT[''.join(sorted([TOLC[y['res']] for y in x]))]) for x in atomlist]


# test code below here
#a = filter_target(parse_pdb_file(open_pdb_file('/home/asriniva/Downloads/3KXU.pdb'), ['ATOM', 'REMARK'], [465]), chain=['A'], name=['CA'])
#b = tessellate(a).vertices
#c = [[a[x] for x in y] for y in b]
#print(c)
#[print(x) for x in potential(c)]
