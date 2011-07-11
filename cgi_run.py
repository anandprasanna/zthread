#!/usr/bin/python
# CGI script for Delaunay Tessellator
# By Arjun Srinivasan
# Jul 19 2010

import cgi, cgitb, urllib2

from PTessellation import PTessellation
from Bio.PDB import PDBParser

def init():
    print "Content-Type: text/html" 
    print
    print "<html><head><title>Output</title></head><body>"


def main():
    cgitb.enable()
    form = cgi.FieldStorage()
    task_name = form["task_name"].value
    task_parse = form["task_parse"].value
    task_tessellate = form["task_tessellate"].value
    task_thread = form["task_thread"].value
    PDB_download = form["PDB_download"].value
    PDB_id = form["PDB_id"].value
    chain = form["chain"].value

    p = PTessellation(PDBParser().get_structure('a', urllib2.urlopen(pdburl+pdbname+".pdb", 'rU')), chain, [srange, erange, irange])
    p.execute()
    if program == "simplex":
        p.print_tessellation()
    elif program == "pprofile":
        p.set_potential(dict([(x[:4], float(x.split()[1])) for x in map(str.rstrip, open(sys.argv[1], 'r').readlines())]))
        p.print_potentials()
    print "</body></html>"

#cgi.test()

if __name__ == '__main__':
    #	cgi.test()
    init()
    main()
