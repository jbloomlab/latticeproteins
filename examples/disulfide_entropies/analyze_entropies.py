"""Analyzes entropies output by disulfide_entropies.py.

Jesse Bloom, 2011."""


import math
import re
import latticeproteins.conformations


def Entropy(a, r, base=20):
    """Returns site entropy at position r (1, 2, 3, ...)"""
    counts = {}
    for (head, seq) in a:
        aa = seq[r - 1].upper()
        if aa in counts:
            counts[aa] += 1
        else:
            counts[aa] = 1
    h = 0.0
    for (aa, n) in counts.iteritems():
        p = n / float(len(a))
        h -= p * math.log(p, base)
    return h


def main():
    """Main body of script."""
    conformations = [
                     'UUURRULURRDDDDDLLUURD',
                    ]
    wo_match = re.compile('^\# Starting with sequence (?P<seq>[A-Z]{22}) without disulfides$')
    w_match = re.compile('^\# Starting with sequence (?P<seq>[A-Z]{22}) with disulfides$')
    for conf in conformations:
        lines = open('%s.log' % conf).readlines()
        initial_seqs = {}
        entropies = {}
        for (match, info) in [(wo_match, 'without'), (w_match, 'with')]:
            for i in range(len(lines)):
                m = match.search(lines[i])
                if m:
                    break
            else:
                raise ValueError("Failed to find sequence.")
            initial_seqs[info] = m.group('seq')
            seqs = []
            for line in lines[i + 1 : ]:
                if line[0] == '#':
                    break
                seqs.append(line.strip())
#            assert len(seqs) == 500
            length = len(seqs[0])
            assert length == 22
            entropies[info] = [Entropy([('head', iseq) for iseq in seqs], r) for r in range(1, length + 1)]
        print "\n---------------\nFor conformation %s" % conf
        print "\nHere is the protein without disulfides:"
        latticeproteins.conformations.PrintConformation(initial_seqs['without'], conf)
        print "\nHere is the protein with disulfides:"
        latticeproteins.conformations.PrintConformation(initial_seqs['with'], conf)
        print "\nRESIDUE\tWITHOUT_DISULFIDES_ENTROPY\tWITH_DISULFIDES_ENTROPY\tENTROPY_DIFFERENCE"
        for r in range(length):
            print "%d\t%.3f\t%.3f\t%.3f" % (r + 1, entropies['without'][r], entropies['with'][r], entropies['with'][r] - entropies['without'][r])


main() # run the script
