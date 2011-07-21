#/usr/bin/env python2
# PDB Parser
# Written by Arjun Srinivasan

import sys
import BeautifulSoup
import itertools
from urllib import urlopen
from re import findall, MULTILINE
from numpy import array
import re
# url for pdb file downloads
PDBFTP = "http://pdbbeta.rcsb.org/pdb/files/"
#url for dssp seq
DSSPFTP = "http://www.pdb.org/pdb/explore/sequenceText.do?structureId="
#specified chain
DSSPCHAIN = "&chainId="
# file extension of pdb files
PDBEXT = ".pdb"

def download_pdb_file(pdbid):
    """ Downloads a pdb file from the Protein Data Bank (PDB)

    Arguments:
    pdbid -- The PDB id of the pdb file.

    Returns:
    File-like url object.
    """
    return urlopen(PDBFTP + pdbid + PDBEXT)

def open_pdb_file(pdbfile):
    """ Opens a locally stored PDB file

    Arguments:
    pdbfile -- The filename of the PDB file.

    Returns:
    File object of the PDB file.
    """
    return open(pdbfile, 'r')

def parse_pdb_file(pdbfile, handle='all', remark_handle='all'):
    """ Parses a PDB file line-by-line, sending individual lines to dedicated parsers.

    Arguments:
    pdbfile -- File object of PDB file.
    handle -- Specific records to handle; if handle="all", all records will be handled.

    Returns:
    Map of PDB file.
    """
    pdbdata = [y for y in map(lambda x: str(x), pdbfile.readlines())]
    pdbfile.close()
    # uses eval to evaluate specific records, using regexes to find said records if the records are in handle.
    return {z: eval("handle_"+z.lower()+"_record("+str(findall(r'^'+z+r'.*\n', ''.join(pdbdata), MULTILINE))+(','+str(remark_handle) if z == 'REMARK' else '')+")") for z in {x[0:6].rstrip() for x in pdbdata} if z in handle or handle=='all'}


def handle_atom_record(records):
    """ Handles atom records.

    Arguments:
    records -- All atom records.

    Returns:
    A list of dictionaries of related information.
    """
    return [{'record' : 'ATOM', 'serial' : int(x[6:11].strip()), 'name' : x[12:16].strip(), 'altloc' : x[16], 'res' : x[17:20], 'chain' : x[21], 'resseq' : int(x[22:26].strip()), 'icode' : x[26], 'coord' : array([float(x[30:38].rstrip()), float(x[38:46].rstrip()), float(x[46:54].rstrip())]), 'occupancy' : float(x[54:60].rstrip()), 'temp' : float(x[60:66].rstrip()), 'element' : x[76:78].lstrip(), 'charge' : x[78:80]} for x in records]

def handle_dbref_record(records):
    """ Handles dbref records.

    Arguments:
    records -- The dbref record.
    """
    return {'record' : 'DBREF', 'id' : records[0][7:11], 'chain' : records[0][12], 'seqbegin' : int(records[0][14:18].strip()), 'insertbegin' : records[0][18], 'seqend' : int(records[0][20:24].strip()), 'insertend' : records[0][24], 'database' : records[0][26:32].strip(), 'dbaccess' : records[0][33:41], 'dbid' : records[0][42:54], 'dbseqbegin' : int(records[0][55:60].strip()) , 'idbnsbegin' : records[0][60], 'dbseqend' : int(records[0][62:67].strip()), 'dbinsend' : records[0][67]}


def handle_remark_record(records, remark_handle):
    """ Handles remark records by passing them on to helper functions for individual remarks.

    Arguments:
    records -- All remark records.

    Returns:
    A dictionary containing information for each type of remark.
    """
    # essentially identical in idea to parse_pdb_file, except it does it for remark numbers.
    return {int(z): eval("handle_remark_"+z+"_record("+str(findall(r'^REMARK'+r'\s*'+z+r'.*\n', ''.join(records), MULTILINE))+")") for z in {x[6:10].strip() for x in records} if int(z) in remark_handle}

def handle_remark_465_record(records):
    """ Handles remark 465, missing residues record.

    Arguments:
    records -- Missing residue records.

    Returns:
    List of missing residues and their numbers (tuple).
    """
    return [{'res' : x[15:18], 'chain' : x[19], 'resseq' : int(x[21:26].strip()), 'icode' : x[26]} for x in records[7:]]
def make_secondary_struc_dic(pdbid,chainid):
   # f=open('output.txt','w')
   # text = ''.join(BeautifulSoup.BeautifulSoup(urlopen(DSSPFTP+pdbid+DSSPCHAIN+chainid)).body(text=True))
    #listq = []
    text2 = ''.join(BeautifulSoup.BeautifulSoup(urlopen(DSSPFTP+pdbid+DSSPCHAIN+chainid)).body(text=True)).replace("&nbsp;","C")
    text2 = re.sub("\d","",text2)# removes all digits
    text2 = re.sub("\n\s*\n*", "\n", text2).strip()# removed multimple newlines replaces with one
    listq = (text2.split("\n"))
    listq.pop()
    listq.pop(0)
    listq = [x for x in listq if listq.index(x) %2 != 0]
    for index, item in enumerate(listq[:]):
        listq[index] = list(item)
    for item in listq:
        y =10
        for x in xrange(10,len(item),y):
            y-=1
            item.pop(x)


    return list(itertools.chain.from_iterable(listq))







"""
    for item in ''.join(BeautifulSoup.BeautifulSoup(urlopen(DSSPFTP+pdbid+DSSPCHAIN+chainid)).body(text=True)).replace("&nbsp;","C").split("\n"):
        print repr(item)
        if((item != '' or item != ' ')and  not (item.isdigit())):
            listq.append(list(item))
"""







# test code below here
#if __name__ == '__main__':
 #   print(parse_pdb_file(open_pdb_file(sys.argv[1]), ['ATOM','REMARK', 'DBREF'], [465]))
