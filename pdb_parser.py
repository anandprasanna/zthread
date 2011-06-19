#/usr/bin/env python3
# PDB Parser
# Written by Arjun Srinivasan

from urllib.request import urlopen
from re import findall, MULTILINE

# url for pdb file downloads
PDBFTP = "http://pdbbeta.rcsb.org/pdb/files/"

# file extension of pdb files
PDBEXT = ".pdb"

def download_pdb_file(pdbid):
	''' Downloads a pdb file from the Protein Data Bank (PDB)
	
	Arguments:
	pdbid -- The PDB id of the pdb file.

	Returns:
	File-like url object.
	'''
	return urlopen(PDBFTP + pdbid + PDBEXT)

def open_pdb_file(pdbfile):
	''' Opens a locally stored PDB file

	Arguments:
	pdbfile -- The filename of the PDB file.

	Returns:
	File object of the PDB file.
	'''
	return open(pdbfile, 'r') 

def parse_pdb_file(pdbfile, handle='all', remark_handle='all'):
	''' Parses a PDB file line-by-line, sending individual lines to dedicated parsers.

	Arguments:
	pdbfile -- File object of PDB file.
	handle -- Specific records to handle; if handle="all", all records will be handled.

	Returns:
	Map of PDB file.
	'''
	pdbdata = pdbfile.readlines()
	pdbfile.close()
	# uses eval to evaluate specific records, using regexes to find said records if the records are in handle.
	return {z: eval("handle_"+z.lower()+"_record("+str(findall(r'^'+z+r'.*\n', ''.join(pdbdata), MULTILINE))+(','+str(remark_handle) if z == 'REMARK' else '')+")") for z in {x[0:6].rstrip() for x in pdbdata} if z in handle or handle=='all'}
	

def handle_atom_record(records):
	''' Handles atom records.

	Arguments:
	records -- All atom records.

	Returns:
	A list of dictionaries of related information.
	'''
	return [{'record' : 'ATOM', 'serial' : int(x[6:11].rstrip()), 'name' : x[12:16].rstrip(), 'altloc' : x[16], 'res' : x[17:20], 'chain' : x[21], 'resseq' : int(x[22:26].rstrip()), 'icode' : x[26], 'coord' : (float(x[30:38].rstrip()), float(x[38:46].rstrip()), float(x[46:54].rstrip())), 'occupancy' : float(x[54:60].rstrip()), 'temp' : float(x[60:66].rstrip()), 'element' : x[76:78].lstrip(), 'charge' : x[78:80]} for x in records]


def handle_remark_record(records, remark_handle):
	'''Handles remark records by passing them on to helper functions for individual remarks.

	Arguments:
	records -- All remark records.

	Returns:
	A dictionary containing information for each type of remark.
	'''
	# essentially identical in idea to parse_pdb_file, except it does it for remark numbers.
	return {int(z): eval("handle_remark_"+z+"_record("+str(findall(r'^REMARK'+r'\s*'+z+r'.*\n', ''.join(records), MULTILINE))+")") for z in {x[6:10].strip() for x in records} if int(z) in remark_handle}

def handle_remark_465_record(records):
	''' Handles remark 465, missing residues record.

	Arguments:
	records -- Missing residue records.

	Returns:
	List of missing residues and their numbers (tuple).
	'''
	return records


# test code below here
print(parse_pdb_file(open_pdb_file('/home/asriniva/Downloads/3KXU.pdb'), ['REMARK'], [465]))
