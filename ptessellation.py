#!/usr/bin/env python2
# Delaunay Tessellator
# Written by Arjun Srinivasan

import sys
from numpy import array,sum as arrsum
from scipy.spatial import Delaunay
from pdb_parser import parse_pdb_file, download_pdb_file, open_pdb_file

def filter_target(atoms, name='all', span='all', chain='all',altloc='A'):
    """ Filters a protein structure based on given information.

    Arguments:
    atoms -- A list of atom records as created by the parser.
    name -- A list of atom types to handle. For example, CA.
    span -- Resseq values that are accepted. Values are entered as arguments to range.
    chain -- Chains to parse.
    altloc -- Specify alternate location to use when faced with a choice.

    FIXME: need to fix find_optimal_span and its usage with filter_target
    """

    return [x for x in atoms if (span == 'all' or x['resseq'] in span) and (name == 'all' or x['name'] in name) and (chain == 'all' or x['chain'] in chain) and (x['altloc']== ' ' or x['altloc'] == altloc) ]


def tessellate(atoms):
    """ Performs a Delaunay tessellation on a given list of atoms.

    Arguments:
    atoms -- List of atoms to triangulate.

    Returns:
    A SciPy Delaunay object (contains info about vertices and neighbors)
    """
    return Delaunay(array([x['coord'] for x in atoms]))


def find_optimal_span(atoms):
    """ Finds the optimal region of the pdb structure, without gaps.

    Arguments:
    atoms -- List of atoms from the pdb structure.

    Returns:
    List of residue numbers of the optimal span.
    """
    # Warnings for missing residues:
#	if 'REMARK' in pdbstruct and 465 in pdbstruct['REMARK']:
#		print("The following residues are missing (residue, chain, resseq, icode): "+str(pdbstruct['REMARK'][465]))
    end_index, max_l = 0,0
    cur_l = 0
    res_seq_numbers = sorted([x['resseq'] for x in atoms])
    for x,y in zip(res_seq_numbers[:-1],res_seq_numbers[1:]):
        if y-x == 1:
            cur_l+=1
        else:
            if max_l < cur_l+1:
                max_l = cur_l+1
                end_index = x
    if max_l < cur_l+2:
        max_l = cur_l+2
        end_index = res_seq_numbers[-1]
    return res_seq_numbers[end_index-max_l:end_index] 




# test code below here
if __name__ == '__main__':
    a = filter_target(parse_pdb_file(open_pdb_file(sys.argv[1]), ['ATOM', 'REMARK'], [465]), chain=['A'], name=['CA'])
    b = tessellate(a).vertices
    c = [[a[x]['res'] for x in y] for y in b]
    d = simplex_potential(c)
#	[print(x['resseq'],x['res'],y) for x,y in zip(a,residue_potential(a, b, d))]

