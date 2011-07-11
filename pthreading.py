#!/usr/bin/env python2
# Simple Protein Threader
# Written by Arjun Srinivasan

import sys
from ptessellation import filter_target
from pprofile import T3LC
from pdb_parser import open_pdb_file, parse_pdb_file

from numpy import array
#from fasta_parser import open_fasta_file, parse_fasta_file

def thread_sequence(mut_seq, orig_seq, orig_seq_begin, seq_numbers, span='all', replace=True):
    """ Threads a protein sequence through a pdb structure. 

    Arguments:
    mut_seq -- Sequence to be threaded.
    orig_seq -- Primary sequence, to be used for replacement.
    orig_seq_begin -- Starting number for PDB numbering.
    span -- Range of residues to use.
    replace -- Flags the option to replace all gaps in sequences with residues from the original template structure.
    """
    return zip(*[(j, T3LC[x] if (x!='-' and x!='X') else T3LC[orig_seq[i]]) for i,x,j in zip(range(len(mut_seq)), mut_seq,seq_numbers) if span=='all' or i+orig_seq_begin in span])


def trim_gaps(mut_seq, orig_seq, seq_numbers):
    """ Trims gaps of the primary sequence.

    Arguments:
    mut_seq -- Sequence to be trimmed. 
    orig_seq -- Primary sequence to be trimmed of gaps.
    seq_numbers -- Numbering for the sequence.
    
    Returns:
    List of sequence numbers and the sequence itself, both trimmed.
    """
    return zip(*[(i,y) for i,x,y in zip(seq_numbers, orig_seq, mut_seq) if x != '-'])

def number_sequence(seq):
    """ Numbers a sequence (as in, doesn't number the gaps).

    Arguments:
    seq -- Sequence to be numbered.

    Returns:
    List of numbers corresponding to each index of the MSA.
    """
    numbers = []
    count = 1
    for x in seq:
        if x != '-':
            numbers+=[count]
            count+=1
        else:
            numbers+=['-']
    return numbers

# main code below here
if __name__ == '__main__':
    a = parse_pdb_file(open_pdb_file(sys.argv[1]), ['ATOM', 'DBREF', 'REMARK'], [465])
    a['ATOM'] = filter_target(a, name=['CA'])
    b = parse_fasta_file(open_fasta_file(sys.argv[2]))
    print(b)
    print(thread_sequence(b[0][1], a))
