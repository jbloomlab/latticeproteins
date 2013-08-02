"""Script for making the topology plot for neato Graphviz program.

Jesse Bloom, 2011."""


def main():
    """Main body of script."""
    
    # input variables
    targets = [ # tuples: outfile, conformation, entropy file (or None), disulfides [(start, end)] numbered 1, ...
                     ('with_disulfides.dot', 'UUURRULURRDDDDDLLUURD', 'UUURRULURRDDDDDLLUURD_entropy_with_disulfides.txt', [(6, 13)]),
                     ('without_disulfides.dot', 'UUURRULURRDDDDDLLUURD', 'UUURRULURRDDDDDLLUURD_entropy_without_disulfides.txt', []),
                     ('targetconformation.dot', 'UUURRULURRDDDDDLLUURD', None, []),
                     ('unfolded1.dot', 'URRDDRRRDLLLLULDDRDRU', None, []),
                     ('unfolded2.dot', 'URDRRDDLDRRDLDLULULUU', None, []),
                     ('unfolded3.dot', 'UURULLDDDDRDDLDLULURU', None, []),
                     ('unfolded4.dot', 'UUUURUURRUURRDDDRDDLL', None, []),
                    ]
    (entropy_max, entropy_min) = (0.54, 0.0) # hardcoded min/max, if None calculated for that protein
    nodeheight = 1 # height of nodes
    bondwidth = 15 # pen width of bonds
    posscale = 100 # scale position values by this much

    for (outfile, conformation, entropyfile, disulfides) in targets:
        # make the graph file
        if entropyfile:
            entropies = dict([(int(x.split()[0]), float(x.split()[1])) for x in open(entropyfile).readlines()])
        else:
            entropies = None

        def EntropyToColor(r):
            if not entropies:
                return 'black'
            else:
                h = entropies[r]
                if entropy_max != None:
                    (max_h, min_h) = (entropy_max, entropy_min)
                else:
                    max_h = max(entropies.values())
                    min_h = min(entropies.values())
                return '"%f, 0.8, 1"' % (1.0 - (max_h - h) / (max_h - min_h) / 3.)

        f = open(outfile, 'w')
        f.write('graph G { layout=neato;\n')
        (x, y) = (0, 0)
        r = 1
        f.write('\tr%d [shape=circle, height=%f, label="", color=%s, style=filled, pin=true, pos="%f,%f"];\n' % (r, nodeheight, EntropyToColor(r), x * posscale, y * posscale))
        for bond in conformation:
            x += {'U':0, 'R':1, 'D':0, 'L':-1}[bond]
            y += {'U':1, 'R':0, 'D':-1, 'L':0}[bond]
            f.write('\tr%d [shape=circle, height=%f, label="", color=%s, style=filled, pin=true, pos="%f,%f"];\n' % (r + 1, nodeheight, EntropyToColor(r + 1), x * posscale, y * posscale))
            f.write('\tr%d -- r%d [penwidth=%f, weight=0];\n' % (r, r + 1, bondwidth))
            r += 1
        for (r1, r2) in disulfides:
            f.write('\tr%d -- r%d [penwidth=%f, weight=0, color=green];' % (r1, r2, bondwidth))
    f.write('}\n')
    f.close()



main() # run the script.
