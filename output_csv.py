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

def write_threading_output(writer, wild_atoms, mut_atoms, wild_potentials, mut_potentials):
    """ Writes output from threading."""
    [writer.writerow(w['resseq'], w['res'], x, y['resseq'], y['res'], z) for w,x,y,z in zip(wild_atoms, wild_potentials, mut_atoms, mut_potentials)]
