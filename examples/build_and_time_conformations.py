"""Script builds the conformations databases, and tests how long it takes to fold proteins.

Written by Jesse Bloom, 2011."""

import time
import os
import latticeproteins.sequences
import latticeproteins.conformations
import latticeproteins.interactions

def main():
    """Main body of script."""
    database_dir = "%s/database" % os.path.split(os.getcwd())[0]    
    ntest = 500 # time by folding this many conformations
    for length in range(10, 25):
        start = time.time()
        print "\nChecking for conformations for protein of length %d, creating if they do not exist." % length
        c = latticeproteins.conformations.Conformations(length, database_dir)
        end = time.time()
        print "Conformations already existed or have been created."
        print "Total time elapsed: %.1f seconds" % (end - start)
        print "Now folding %d proteins to calculate the time to fold." % ntest
        start = time.time()
        for i in range(ntest):
            s = latticeproteins.sequences.RandomSequence(length)
            (dg, conf, numcontacts) = c.FoldSequence(s, 1.0)
        end = time.time()
        print "Took a total of %.3f seconds to fold %d sequences, or %.3f seconds per sequence." % (end - start, ntest, (end - start) / float(ntest))
        print "The total number of contact sets is %d." % c.NumContactSets()



main() # run the script.
