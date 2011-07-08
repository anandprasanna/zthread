#!/usr/bin/env python2
# CSV Output Program
# Written by Arjun Srinivasan

import csv

def open_csv_file(filename):
    """ Open CSV file for writing.

    Arguments:
    filename -- CSV file.

    Returns:
    CSV writer object.
    """
    return csv.writer(open(filename, 'w'), delimiter=' ')

def write_simplex_potential():
    """ Writes information about individual simplices."""
    pass

def write_residue_potential():
    """ Writes information about residue potentials. """
    pass

def write_threading_output(writer, wild_numbers, wild_residues, wild_potentials, mut_numbers, mut_residues, mut_potentials):
    """ Writes output from threading."""
    for u,v,w,x,y,z in zip(wild_numbers, wild_residues, wild_potentials, mut_numbers, mut_residues, mut_potentials): writer.writerow([u,v,w,x,y,z])
