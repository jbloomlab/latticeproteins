"""Script for examining conformations visually.

Written by Jesse Bloom, 2011."""

import os
import random
import latticeproteins.sequences
import latticeproteins.conformations
import latticeproteins.interactions

def main():
    """Main body of script."""
    database_dir = "%s/database" % os.path.split(os.getcwd())[0]    
    length = 22
    c = latticeproteins.conformations.Conformations(length, database_dir)
    seq = ''.join(['A' for x in range(length)])
    uniqueconfs = []
    for i in range(c.MaxContacts() - 4, 1, -1):
        uniqueconfs += c.UniqueConformations(i)
    random.shuffle(uniqueconfs)
    print "Randomly stepping through the %d unique conformations one by one." % len(uniqueconfs)
    for conf in uniqueconfs:
        x = raw_input("Press ENTER for next conformation.")
        print "\nConformation %s" % conf
        latticeproteins.conformations.PrintConformation(seq, conf)



main() # run the script.
