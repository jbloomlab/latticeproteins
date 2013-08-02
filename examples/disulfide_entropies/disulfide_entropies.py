"""Examines how disulfide bonds influenza site entropies.

Written by Jesse Bloom, 2011."""

import sys
import random
import latticeproteins.sequences
import latticeproteins.conformations
import latticeproteins.interactions
import latticeproteins.fitness
import latticeproteins.evolution


def main():
    """Main body of script."""
    random.seed(1) # seed random number generator for predictable behavior
    database_dir = "/Users/jbloom/latticeproteins-0.1/database" 
    targets = [
               ('UUURRULURRDDDDDLLUURD', 6, 13),
              ]
    length = len(targets[0][0]) + 1 # protein length
    print "Loading conformations for length %d" % length
    confs = latticeproteins.conformations.Conformations(length, database_dir, interaction_energies=latticeproteins.interactions.miyazawa_jernigan_disulfide) # holds conformations
    maxsteps = 500
    f_cutoff = 0.1
    for (targetconf, disulfide_start, disulfide_end) in targets:
        assert length == len(targetconf) + 1
        print '\nFor target conformation %s.' % targetconf
        print "Creating fitness evaluator."
        f_threshold = latticeproteins.fitness.Fitness(1.0, confs, 0.0, targetconf) # sequences that fold to targetconf with dG <= 0 have fitness 1, otherwise fitness of 0
        f_linear = latticeproteins.fitness.Fitness(1.0, confs, 'negstability', targetconf) # fitness is negative of stability
        print "Performing an adaptive walk to find a folding sequence." 
        seq = ''.join(latticeproteins.sequences.RandomSequence(length)).replace('C', 'A')
        print "Beginning with the initial sequence of %s." % seq
        (seq, f, n) = latticeproteins.evolution.AdaptiveWalk(seq, f_linear, maxsteps, f_cutoff, exclude_mutations=dict([(r, 'C') for r in range(length)]))
        if f >= f_cutoff:
            print "Found the sequence %s, with stability -%.2f, after %d steps." % (seq, f, n)
        else:
            print "FAILED TO FIND SEQUENCE WITH TARGET STABILITY OF -%.2f AFTER %d STEPS." % (f_cutoff, maxsteps)
            continue
        print "Here is the folded sequence:" 
        latticeproteins.conformations.PrintConformation(seq, targetconf)
        assert 'C' not in seq
        print "\nSubstituting cysteines at base of the loop."
        disulfide_seq = list(seq)
        disulfide_seq[disulfide_start - 1] = 'C'
        disulfide_seq[disulfide_end - 1] = 'C'
        disulfide_seq = ''.join(disulfide_seq)
        latticeproteins.conformations.PrintConformation(disulfide_seq, targetconf)
        print "With these cysteines, the new stability is -%.2f" % f_linear.Fitness(disulfide_seq)
        if f_linear.Fitness(disulfide_seq) < f_cutoff:
            print "Performing an adaptive walk to increase stability to at least %.3f" % f_cutoff
            (disulfide_seq, f, n) = latticeproteins.evolution.AdaptiveWalk(disulfide_seq, f_linear, maxsteps, f_cutoff, exclude_mutations=dict([(r, 'C') for r in range(length) if r not in [disulfide_start - 1, disulfide_end - 1]]))
            print "Here is the folded sequence:"
            latticeproteins.conformations.PrintConformation(disulfide_seq, targetconf)
            if f < f_cutoff:
                print "FAILED: DISULFIDE SEQ DOES NOT STABLY FOLD."
                continue
        assert disulfide_seq.count('C') == 2 and disulfide_seq[disulfide_start - 1] == 'C' and disulfide_seq[disulfide_end - 1] == 'C'
        print "Here is the final sequence with disulfides; stability is -%.2f" % f
        latticeproteins.conformations.PrintConformation(disulfide_seq, targetconf)
        # Now begin the evolution
        print "Now beginning evolution."
        logfile = "%s.log" % targetconf
        f = open(logfile, 'w')
        f.write("# Sequences target conformation %s\n" % targetconf)
        for (iseq, infostring) in [(seq, 'without'), (disulfide_seq, 'with')]:
            f.write("# Starting with sequence %s %s disulfides\n" % (iseq, infostring)) 
            for ireplicate in range(500):
                print "Replicate %d %s disulfides." % (ireplicate, infostring)
                (finalpop, finalstabs, finalsubs) = latticeproteins.evolution.NeutralEvolution(iseq, f_threshold, 0.0, popsize=10, mutrate=0.1 / length, numsteps=100, allproteinsfile='allproteinsfile.temp')
                (seqtosave, stabtosave) = random.choice([(x, y) for (x, y) in zip(finalpop, finalstabs) if y < 0.0 and y != None])
                f.write("%s\n" % seqtosave)
                f.flush()
                print "Saving sequence %s with stability %.3f" % (seqtosave, stabtosave)
                sys.stdout.flush()
        f.close()



main() # run the script.
