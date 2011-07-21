#!/usr/bin/env python2
# Delaunay Tessellator
# Written by Arjun Srinivasan

import sys
from numpy import array,sum as arrsum
from scipy.spatial import Delaunay
from sets import Set
from pdb_parser import parse_pdb_file, download_pdb_file, open_pdb_file
import math
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

    return [x for x in atoms if (span == 'all' or x['resseq'] in span) and (name == 'all' or x['name'] in name) and (chain == 'all' or x['chain'] in chain) and (x['altloc']== ' ' or x['altloc'] == altloc) and x['icode'] == ' ' ]


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
                end_index = y
            cur_l = 0
    if max_l < cur_l+1:
        max_l = cur_l+1
        end_index = len(res_seq_numbers)
    return res_seq_numbers[end_index-max_l:end_index]




def make_matrix(pdbid,chainid):
    matrixlist = []
    makelist = make_matrixlist(pdbid,chainid)
    proxlist = makelist[0]
    simplexlist =  makelist[1]
    for res in proxlist:






















def make_prematrixlist(pdbid,chainid): #makes a matrix given the pdb id and cahin id for now after removing simplicies i do not update data ,i assume that its very unlikly that a residue will have no simplcies with lengths less than 10 a
    pdbstruct = (parse_pdb_file(download_pdb_file(pdbid),['ATOM'],[] ))
    a = filter_target(pdbstruct['ATOM'], name = ['CA'], chain = [chainid])
    a = filter_target(a, name = ['CA'], span=find_optimal_span(a))
    simplexlist = [[ (a[x]['resseq'], a[x]['res'], a[x]['coord'] ) for x in y] for y in tessellate(a).vertices]#gives grouping in terms of four(one simplex)
    data = ([(x['resseq'],x['res']) for x in a])
    simplexlist = filter_simplicies(10,simplexlist)#integer sets cutoff
    return (neighborlist(data,simplexlist),simplexlist)


def neighborlist(reslist,simplexlist):#returns list of neighbors(in set form for each residue in terms of indicies in simnplexlist
    globallist = [None]*(reslist[-1][0]+10)
    for item in reslist:
        globallist[item[0]] = addsecshell(getshells(item[0],simplexlist),simplexlist)
    return globallist

def addsecshell(shell1,simplexlist):#adds second shell thorugh update
    secshell= shell1.copy()
    for simplex in shell1.copy():
        sechell = secshell.update(getshells(simplex,simplexlist))
    return secshell



def getshells(resnum,simplexlist):#returns set
    shellist = Set() #empty set
    for index, item2 in enumerate(simplexlist):
        if(item2[0][0]== resnum or item2[1][0]== resnum or item2[2][0]== resnum or item2[3][0]== resnum):
            shellist.add(index)

    #for secshell in shellist[:]:
    return shellist





def filter_simplicies(cutoff,simplexlist):
    for index,item in enumerate(simplexlist):#checks every simplex for edge lengths greater than cutoff
        if not(calc_dist(item[0][2],item[1][2],cutoff) and calc_dist(item[0][2],item[2][2],cutoff) and calc_dist(item[0][2],item[3][2],cutoff)
            and calc_dist(item[1][2],item[2][2],cutoff) and calc_dist(item[1][2],item[3][2],cutoff) and calc_dist(item[2][2],item[3][2],cutoff)):
            simplexlist.pop(index)

    return simplexlist



def calc_dist(coord1,coord2,cutoff):
    return True if math.sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2 + (coord1[2]-coord2[2])**2)<cutoff else False

#def transform(coordlist):
   # if c



# test code below here
#if __name__ == '__main__':
    a = filter_target(parse_pdb_file(open_pdb_file(sys.argv[1]), ['ATOM', 'REMARK'], [465]), chain=['A'], name=['CA'])
    b = tessellate(a).vertices
    c = [[a[x]['res'] for x in y] for y in b]
    d = simplex_potential(c)
#	[print(x['resseq'],x['res'],y) for x,y in zip(a,residue_potential(a, b, d))]

