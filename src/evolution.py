#!/usr/bin/python
# Begin evolution.py
#---------------------------------------------------------------------------
"""Module for lattice protein evolution algorithms.

Written by Jesse Bloom, 2004."""
#---------------------------------------------------------------------------
import math, conformations, sequences, sys, random, stats
#----------------------------------------------------------------------
class EvolutionError(Exception):
    """Error in lattice protein evolution."""
#----------------------------------------------------------------------
def NeutralEvolution(initialpop, stability, stabilitycutoff, popsize, mutrate, numsteps, allproteinsfile):
    """Performs neutral evolution of lattice proteins.

    Call is: '(finalpop, finalstabs, finalsubs) = NeutralEvolution(initialpop, stability,
        stabilitycutoff, popsize, mutrate, numsteps,
        allproteinsfile)'
    This method neutrally evolves a population of proteins.  
        All sequences either have fitness of one (if 
        'stability[seq] <= stabilitycutoff') or zero (otherwise).
        At each generation, we choose 'popsize' offspring with
        equal probability from all fitness one sequences, and
        use them to replace the old population.  All sequences
        are then mutated with a per site mutation rate of 'mutrate'.
    'initialpop' specifies the initial population.  It can be:
	* a single sequence, in which case the initial population is clonal
    	    copies of this sequence.  If it is a single sequence, it
	    must be a string (not a list)!!!
	* an integer, in which case the population is randomly generated
	    sequences of the length specified by the integer.
	* a list of sequences of length equal to 'popsize', giving
            the initial population.
    'stability' is the object used to evaluate the stability of 
	protein sequences.  It must have a method 'Stability(seq)' 
        that returns the stability of 'seq', or 'None' if
        that sequence is unfolded.
    'popsize' is an integer > 0 specifying the population size.
    'mutrate' is the per residue per generation mutation rate.  
	Each residue is randomly mutated with probability 'mutrate'
	in each generation.  It is a number > 0 and < 1.0.
    'numsteps' is the number of steps for which the evolutionary trajectory
	proceeds.  It is an integer > 0.
    'allproteinsfile' is a string giving the name of a file to which
        we write information on all proteins in the population at
        each generation.
    The returned variable is the 3-tuple '(finalpop, finalstabs, 
        finalsubs)'.  'finalpop' is a list of all of the 
        sequences in the final population, 'finalstabs' is a list
        of the stabilities of all of these sequences, and
        'finalsubs' is a list of the number of substitutions
        of all of these sequences.  If the sequence is unfolded,
        the listed stability may be 'None.'"""
    # do some error checking on the input variables
    if not isinstance(stabilitycutoff, (int, float)):
        raise EvolutionError, "Invalid stability cutoff of %f." % stabilitycutoff
    if not (isinstance(popsize, int) and popsize > 0):
        raise EvolutionError, "Invalid 'popsize' of %r." % popsize
    if not (isinstance(numsteps, int) and numsteps > 0):
        raise EvolutionError, "Invalid 'numsteps' of %r." % numsteps
    if not 0.0 < mutrate < 1.0:
        raise EvolutionError, "Invalid 'mutrate' of %r." % mutrate
    # set the initial population
    if isinstance(initialpop, str): # single sequence, make clonal copies
        initialpop = [initialpop for i in range(popsize)]
    elif isinstance(initialpop, int): # choose random sequences
        initialpop = [sequences.RandomSequence(initialpop) for i in range(popsize)]
    elif isinstance(initialpop, list): # random sequences from this list
        if len(initialpop) != popsize:
            raise EvolutionError, "'initialpop' is of invalid size."
    else:
        raise EvolutionError, "Invalid 'initialpop' of %r." % initialpop
    # print some information
    file = open(allproteinsfile, 'w')
    file.write("Performing a neutral evolutionary run of %d steps.\n" % numsteps)
    file.write("Stability cutoff is %f." % stabilitycutoff)
    file.write("The population size is %d.\n" % popsize)
    file.write("The per site per generation mutation rate is %f.\n" % mutrate)
    file.write("In column headings, 'dG' is free energy of folding, 'F' is\n")
    file.write("the fitness, 'S' is the number of substitutions from ancestor,\n")
    file.write("'N' is the number of these sequences in the population,\n")
    file.write("ID is the ID number of this sequence, and Parent_ID is\n")
    file.write("the ID of the sequence's parent, or 0 if this sequence\n")
    file.write("was in the original population.\n")
    file.write("STEP\tdG\tF\tS\tSequence\tN\tID\tParent_ID\n")
    # population is composed of 5-tuples: 
    #  (sequence, dGf, numsubstitutions, id, parent_id)
    population_dict = {}
    for seq in initialpop:
        seq = ''.join(seq)
        try:
            population_dict[seq] += 1
        except KeyError:
            population_dict[seq] = 1
    idcounter = 1
    population = []
    for (seq, count) in population_dict.iteritems():
        population += [(seq, stability.Stability(seq), 0, idcounter, 0)] * count
        idcounter += 1
    del population_dict
    istep = 0
    # write out the initial population
    prot_dict = {}
    for prot in population:
        try:
            prot_dict[prot] += 1
        except KeyError:
            prot_dict[prot] = 1
    # write info to file
    for ((seq, dGf, numsubs, id, parent_id), n) in prot_dict.iteritems():
        if dGf <= stabilitycutoff and dGf != None:
            ifit = 1.0
        else:
            ifit = 0.0
        if dGf == None:
            file.write("%d\tUNFOLDED\t%.4f\t%d\t%s\t%d\t%d\t%d\n" % (istep, ifit, numsubs, seq, n, id, parent_id))
        else:
            file.write("%d\t%.4f\t%.4f\t%d\t%s\t%d\t%d\t%d\n" % (istep, dGf, ifit, numsubs, seq, n, id, parent_id))
    # begin the evolutionary trajectory
    while True:
        istep += 1
        # mutate all of the sequences in the population
        # also store sequences in 'prot_dict' for writing to file
        # 'prot_dict' is keyed by proteins, values are number in pop.
        prot_dict = {}
        for iprot in range(len(population)):
            prot = population[iprot]
            seq = prot[0]
            mutseq = sequences.MutateSequence(seq, mutrate)
            if mutseq != seq: # protein mutated
                mutseq = ''.join(mutseq)
                mutprot = (mutseq, stability.Stability(mutseq), prot[2] + 1, idcounter, prot[3])
                idcounter += 1
                population[iprot] = mutprot
                try:
                    prot_dict[mutprot] += 1
                except KeyError:
                    prot_dict[mutprot] = 1
            else: # protein not mutated
                prot = (prot[0], prot[1], prot[2], prot[3], prot[3])
                population[iprot] = prot
                try:
                    prot_dict[prot] += 1
                except KeyError:
                    prot_dict[prot] = 1
        # write info to file
        atleastonealive = False
        for ((seq, dGf, numsubs, id, parent_id), n) in prot_dict.iteritems():
            if dGf <= stabilitycutoff and dGf != None:
                ifit = 1.0
                atleastonealive = True
            else:
                ifit = 0.0
            if dGf == None:
                file.write("%d\tUNFOLDED\t%.1f\t%d\t%s\t%d\t%d\t%d\n" % (istep, ifit, numsubs, seq, n, id, parent_id))
            else:
                file.write("%d\t%.4f\t%.1f\t%d\t%s\t%d\t%d\t%d\n" % (istep, dGf, ifit, numsubs, seq, n, id, parent_id))
        if not atleastonealive:
            raise EvolutionError, "No living sequences in population."
        if istep >= numsteps:
            break
        # replace current population with offspring
        newpopulation = []
        while len(newpopulation) < popsize:
            prot = random.choice(population)
            if prot[1] <= stabilitycutoff and prot[1] != None:
                newpopulation.append(prot)
        population = newpopulation
    file.close()
    # make a list of the final population and return
    #  (sequence, dGf, numsubstitutions, id, parent_id)
    finalpop = [prot[0] for prot in population]
    finalstabs = [prot[1] for prot in population]
    finalsubs = [prot[2] for prot in population]
    return (finalpop, finalstabs, finalsubs)
#----------------------------------------------------------------------
def Evolution(initialpop, fitness, popsize, mutrate, numsteps, file = sys.stdout, targetfunction = None, allproteinsfile = None, dGincrement = 0.02, dGlist = None):
    """Performs evolution of lattice proteins.

    Call is: 'finalpop = Evolution(fitness, popsize, mutrate, 
	numsteps, [file = sys.stdout, targetfunction = None, 
	allproteinsfile = None, dGincrement = 0.02, dGlist = {}])'
    This method evolves a population of proteins.  At each generation,
	we choose 'popsize' offspring.  These offspring 
	are identical copies of existing members of the population, with
	the probability of a copy being from a member of the population
	being proportional to the fitness of that member of the 
	population.  These new sequences then replace the 
	existing sequences in the population.  All
	members of the population are then mutated with a per site
	mutation rate of 'mutrate'.  A few points: the probability that
	an offspring is from some parent is proportional to that parent's
	fitness if the fitness is > 0 -- if the fitness is <= zero, then
	the probability that an offspring is from some parent is zero.
	However, if all members of the population have fitnesses <= 0,
	then the offspring are randomly created from these parents
	with equal probability.  Both the selection of parents and 
	the mutation of sequences are stochastic.
    'initialpop' specifies the initial population.  It can be:
	* a single sequence, in which case the initial population is clonal
    	    copies of this sequence.  If it is a single sequence, it
	    must be a string (not a list)!!!
	* an integer, in which case the population is randomly generated
	    sequences of the length specified by the integer.
	* a list of sequences.  If the length of the list is equal to 
	    'popsize', then this list is the initial population.  If
	    the length of this list is greater than 'popsize', then
	    'popsize' sequences from the list are randomly chosen
	    for the initial population.  If the length of the list
	    is less than 'popsize', then additional sequences are 
	    created as identical copies of randomly selected members of
	    the list.
    'fitness' is the object used to evaluate the fitness of the protein
	sequences.  It must have a method 'Fitness(seq)' such that
	'fitness.Fitness(seq)' returns a number representing the fitness
	of protein sequence 'seq', where 'seq' is a member of the population.
	as 'initialseq'.  If 'file' is not 'None', it must also have a 
	method 'Info(file)' that prints information about the fitness
	evaluation to 'file'.  Fitnesses are numbers, larger numbers
	imply higher fitness.  Note that any fitnesses less than zero
	are set to zero for the purpose of specifying reproduction rate.
    'popsize' is an integer > 0 specifying the population size.
    'mutrate' is the per residue per generation mutation rate.  Each residue
	in each sequence is randomly mutated with probability 'mutrate'
	in each generation.  It is a number > 0 and < 1.0.
    'numsteps' is the number of steps for which the evolutionary trajectory
	proceeds.  It is an integer > 0.
    'file' is an optional argument specifying where the output information
	for the walk is printed.  The default value is standard output
	('file = sys.stdout').  If 'file' is set to an open writable
	file-like object, then output is printed to that object instead.
	If 'file' is set to 'None', then no output is printed.
    'targetfunction' is an optional argument that tells us to stop when 
	the evolutionary run reaches some target, rather than going for
	'numsteps' steps.  By default, 'targetfunction' is 'None'.  To
	set this option, set 'targetfunction' to the 3-tuple
	'(function, targetvalue, range)'.  'function' is a method that
	takes as input a single protein sequence and returns a number.
	'targetvalue' is a number, and 'range' is a positive number.
	As soon as the evolutionary run generates a sequence 'seq'
	such that 'function(seq)' is within 'functionrange' of 'targetvalue',
	the evolutionary run stops and 'seq' is returned as 'finalpop'.
	If the target is not reached in 'numsteps' steps, an exception
	is raised.
    'allproteinsfile' is an optional argument that specifies that we write
	data giving the details of all proteins in the population at each
	step of evolution.  The default value is 'None', meaning nothing
	is written.  If it is set to a non-'None', it should be set to
	a writable file-like object.  The data is then written to
	this file, along with an explanatory header.
    'dGincrement' and 'dGlist' allow the program to keep track of 
        the stabilities of all of the proteins generated during the
        course of the evolutionary run.  This is done by passing
        a variable representing an empty list as 'dGlist'.
        This variable is not explicitly returned by the function,
        but when the function is completed 'dGlist' will contain
        information about the stabilities of the sequences generated.
        This is essentially histogram-style information, with a bin
        size of 'dGincrement'.  'dGlist[0]' is the number of sequences
        generated having stabilities dG such that dG > 0. For dG < 0,
        'dGlist[i]' is the number of sequences generated having 
        stabilities dG such that:
            dGincrement * (i - 1) <= dG < dGincrement * i
        Note that the only valid value for 'dGlist' is either the
        default value of 'None' or a variable representing an empty
        list, and the 'dGincrement' must be greater than zero.
    The returned variable is the list 'finalpop' which gives all of the
	sequences in the final population.  These sequences are sorted
	from highest fitness to lowest fitness.  However, if 
	'targetfunction' is set to a non-'None' value (see above) only
	a single sequence is returned."""
    # do some error checking on the input variables
    if not (isinstance(popsize, int) and popsize > 0):
	raise EvolutionError, "Invalid 'popsize' of %r." % popsize
    if not (isinstance(numsteps, int) and numsteps > 0):
	raise EvolutionError, "Invalid 'numsteps' of %r." % numsteps
    if not 0.0 < mutrate < 1.0:
	raise EvolutionError, "Invalid 'mutrate' of %r." % mutrate
    if targetfunction == None:
	pass
    elif isinstance(targetfunction, tuple) and len(targetfunction) == 3:
	(function, targetvalue, functionrange) = targetfunction
	if not callable(function):
	    raise EvolutionError, "'targetfunction[0]' is not callable."
	if not isinstance(targetvalue, (int, float)):
	    raise EvolutionError, "'targetfunction[1]' is not a number."
	if not isinstance(functionrange, (int, float)) and functionrange >= 0.0:
	    raise EvolutionError, "'targetfunction[2] is not a positive number."
    else:
	raise EvolutionError, "Invalid 'targetfunction' of %r." % targetfunction
    if not (dGlist == [] or dGlist == None):
        raise EvolutionError, "Invalid 'dGlist' of %r." % dGlist
    if dGlist != None:
        dGlist.append(0)
    dGincrement = float(dGincrement)
    if not (isinstance(dGincrement, float) and dGincrement > 0):
        raise EvolutionError, "Invalid 'dGincrement' of %r." % dGincrement
    # set the initial population
    if isinstance(initialpop, str): # single sequence, make clonal copies
	population = [initialpop for i in range(popsize)]
    elif isinstance(initialpop, int): # choose random sequences
	population = [sequences.RandomSequence(initialpop) for i in range(popsize)]
    elif isinstance(initialpop, list): # random sequences from this list
	if len(initialpop) == popsize:
	    population = list(initialpop)
	elif len(initialpop) > popsize:
	    population = list(initialpop)
	    random.shuffle(population)
	    population = population[0 : popsize]
	else:
	    population = list(initialpop)
	    while len(population) < popsize:
		population.append(random.choice(population))
    else:
	raise EvolutionError, "Invalid 'initialpop' of %r." % initialpop
    # print some information
    if file != None:
	file.write("Performing an evolutionary run of %d steps.\n" % numsteps)
	file.write("The population size is %d.\n" % popsize)
	file.write("The per site per generation mutation rate is %f.\n" % mutrate)
	if targetfunction:
	    file.write("Evolving to a target.  Evolution continues until the value\n")
	    file.write("returned by %r is within %f of %f." % (function, functionrange, targetvalue))
	fitness.Info(file)
	file.write("In column headings, 'F' is fitness, 'dG' is free energy of\n")
	file.write("folding, and 'S' is the number of substitutions from ancestor.\n")
	file.write("'<X>' is the mean of X, 'MED_X' is the median of X,\n")
	file.write("'BEST_X' is the value of X for the most fit sequence,\n") 
	file.write("'MOST_X' is the value of X for the most abundant sequence.\n")
	file.write("STEP\tBEST_F\t<F>\tMED_F\tMOST_F\tBEST_dG\t<dG>\tMED_dG\tMOST_dG\tBEST_S\t<S>\tMED_S\tMOST_S\tBEST_SEQ\n")
    if allproteinsfile != None:
	allproteinsfile.write("Performing an evolutionary run of %d steps.\n" % numsteps)
	allproteinsfile.write("The population size is %d.\n" % popsize)
	allproteinsfile.write("The per site per generation mutation rate is %f.\n" % mutrate)
	if targetfunction:
	    allproteinsfile.write("Evolving to a target.  Evolution continues until the value\n")
	    allproteinsfile.write("returned by %r is within %f of %f." % (function, functionrange, targetvalue))
	fitness.Info(allproteinsfile)
	allproteinsfile.write("In column headings, 'F' is fitness, 'dG' is free energy of\n")
	allproteinsfile.write("folding, 'S' is the number of substitutions from ancestor,\n")
	allproteinsfile.write("and 'N' is the number of these sequences in the population.\n")
        allproteinsfile.write("ID is the ID number of this sequence, and Parent_ID is\n")
        allproteinsfile.write("the ID of the sequence's parent, or 0 if this sequence\n")
        allproteinsfile.write("was in the original population.\n")
        allproteinsfile.write("STEP\tdG\tF\tS\tSequence\tN\tID\tParent_ID\n")
    # keep track of substitutions of each sequence from its ancestor
    substitutions = [0] * popsize
    # assign IDs to sequences and assign parent IDs
    parent_ids = [0] * popsize
    ids = [None] * popsize
    idcounter = 1
    seq_index_dict = {}
    for i in range(popsize):
        seq = ''.join(population[i])
        try:
            seq_index_dict[seq].append(i)
        except KeyError:
            seq_index_dict[seq] = [i]
    for l in seq_index_dict.itervalues():
        for i in l:
            ids[i] = idcounter
        idcounter += 1
    assert not ids.count(None)
    # begin the evolutionary trajectory
    istep = 0
    while True:
	# shuffle the population and substitutions lists
	combined_list = zip(substitutions, population, ids, parent_ids)
	random.shuffle(combined_list)
	substitutions = [tup[0] for tup in combined_list]
	population = [tup[1] for tup in combined_list]
        ids = [tup[2] for tup in combined_list]
        parent_ids = [tup[3] for tup in combined_list]
	# evaluate the fitnesses
	fitnesses = [max(0.0, fitness.Fitness(seq)) for seq in population]
	# evaluate the stabilities 
	stabilities = [fitness.Stability(seq) for seq in population]
	# compute the best fitness and the "total fitness"
	bestfit = fitnesses[0]
	bestseq = population[0]
	beststab = stabilities[0]
	bestsubs = substitutions[0]
	totalfitness = 0.0
	for i in range(len(fitnesses)):
	    ifit = fitnesses[i]
	    totalfitness += ifit
	    if ifit > bestfit:
		bestfit = ifit
		bestseq = population[i]
		beststab = stabilities[i]
		bestsubs = substitutions[i]
            if dGlist:
                istab = stabilities[i]
                if istab > 0:
                    dGlist[0] += 1
                else:
                    dGlist_index = int(-istab / dGincrement) + 1
                    try:
                        dGlist[dGlist_index] += 1
                    except IndexError:
                        dGlist += [0] * (dGlist_index + 1 - len(dGlist))
                        dGlist[dGlist_index] += 1
	# write information to file
	if file != None:
	    # find the most abundant sequence
	    most_seq = sequences.MostAbundant(population)
	    most_index = population.index(most_seq)
	    most_fitness = fitnesses[most_index]
	    most_stability = stabilities[most_index]
	    most_substitutions = substitutions[most_index]
	    # write the information to the file
    	    file.write("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n" % (istep, bestfit, stats.Mean(fitnesses), stats.Median(fitnesses), most_fitness, beststab, stats.Mean(stabilities), stats.Median(stabilities), most_stability, bestsubs, stats.Mean(substitutions), stats.Median(substitutions), most_substitutions, ''.join(bestseq)))
	if allproteinsfile != None:
	    prot_dict = {}
	    for i in range(len(fitnesses)):
		try:
    		    prot_dict[(istep, stabilities[i], fitnesses[i], substitutions[i], ''.join(population[i]), ids[i], parent_ids[i])] += 1
		except KeyError:
    		    prot_dict[(istep, stabilities[i], fitnesses[i], substitutions[i], ''.join(population[i]), ids[i], parent_ids[i])] = 1
	    for (tup, n_in_pop) in prot_dict.iteritems():
	      	allproteinsfile.write("%d\t%.4f\t%.4f\t%.2f\t%s\t%d\t%d\t%d\n" % (tup[0], tup[1], tup[2], tup[3], tup[4], n_in_pop, tup[5], tup[6]))
	# see if we have met the targetfunction
	if targetfunction:
	    for seq in population:
		if math.fabs(function(seq) - targetvalue) <= functionrange:
		    return ''.join(seq)
	# see if we are done
	if istep >= numsteps:
	    if targetfunction:
		raise EvolutionError, "Failed to reach the target function."
	    break
	# choose the offspring 
	cumulative_fitnesses = []
	cumulative_total = 0.0
	for f in fitnesses:
	    try:
	    	x = f / totalfitness
	    except ZeroDivisionError:
		x = 1.0 / popsize
	    cumulative_total += x
	    cumulative_fitnesses.append(cumulative_total)
	newsequences = []
	newsubstitutions = []
        new_ids = []
        new_parent_ids = []
	for inew in range(popsize):
	    x = random.random() # random number from zero to 1
	    if x <= cumulative_fitnesses[0]:
		iseq = 0
	    else:
    		for iseq in range(1, len(cumulative_fitnesses)):
    		    if cumulative_fitnesses[iseq] >= x >= cumulative_fitnesses[iseq - 1]:
    			break
    		else:
    		    if not -1e-5 <= cumulative_fitnesses[iseq] - x <= 1e-5:
    			raise EvolutionError, "Error choosing offspring: %r, %r, %r." % (x, cumulative_fitnesses[iseq], iseq)
	    newsequences.append(population[iseq])
	    newsubstitutions.append(substitutions[iseq])
            new_ids.append(ids[iseq])
            new_parent_ids.append(parent_ids[iseq])
	# replace current population with offspring
	population = newsequences
	substitutions = newsubstitutions
        ids = new_ids
        parent_ids = new_parent_ids
	# mutate all sequences
	for i in range(len(population)):
	    oldseq = population[i]
            newseq = ''.join(sequences.MutateSequence(oldseq, mutrate))
            if oldseq != newseq:
                population[i] = newseq 
                substitutions[i] += sequences.HammingDistance(oldseq, newseq)
                parent_ids[i] = ids[i]
                ids[i] = idcounter
                idcounter += 1
	istep += 1
    # sort the population from highest to lowest fitness
    decorated_list = []
    for i in range(len(population)):
	decorated_list.append((fitnesses[i], population[i]))
    decorated_list.sort()
    decorated_list.reverse()
    population = []
    for x in decorated_list:
	population.append(x[1])
    # return the final population
    return population
#---------------------------------------------------------------------------
def AdaptiveWalk(initialseq, fitness, maxsteps, fitness_cutoff=None, steepest_ascent=False, file=sys.stdout, exclude_mutations=None):
    """Performs an adaptive walk or a steepest ascent walk.

    Call is: '(seq, f, n) = AdaptiveWalk(initialseq, fitness, maxsteps,
	[fitness_cutoff = None, steepest_ascent = False, file = sys.stdout])'
    This method performs an adaptive walk by creating single point
	mutants of a sequence and selecting more fit ones, and then
	repeating this procedure.  The walk 
	can proceed as either an adaptive walk, where the first 
	sequence more fit than wild-type from the previous generation
	is selected, or as a steepest ascent walk, where the point
	mutant that leads to the largest improvement over the wild-type
	is selected.  In either case, the walk terminates when no mutants
	have a higher sequence than wild-type or when the fitness reaches
	some specified cutoff.
    'initialseq' is a string or list specifying the protein sequence
        from which the adaptive walk begins.
    'fitness' is the object used to evaluate the fitness of the protein
        sequences.  It must have a method 'Fitness(seq)' such that
        'fitness.Fitness(seq)' returns a number representing the fitness
        of protein sequence 'seq', where 'seq' is of the same length
        as 'initialseq'.  If 'file' is not 'None', it must also have a 
        method 'Info(file)' that prints information about the fitness
        evaluation to 'file'.  Fitnesses are numbers, larger numbers
        imply higher fitness.
    'maxsteps' is the maximum number of steps in the adaptive walk.
        The walk proceeds either for 'maxsteps' steps, until the next
        step fails to locate a sequence with higher fitness, or until
        the criterium set by 'fitness_cutoff' is achieved.
    'fitness_cutoff' is an optional argument specifying a fitness cutoff.
        The default value of this argument is 'None', if it is set to
        a number, then the walk terminates when the sequence achieves
        a fitness >= 'fitness_cutoff'.
    'steepest_ascent' is an optional argument that specifies that we
        perform a steepest ascent walk rather than an adaptive walk.
        The default is 'False'; if it set to 'True' then we perform
        a steepest ascent walk.
    'file' is an optional argument specifying where the output information
        for the walk is printed.  The default value is standard output
        ('file = sys.stdout').  If 'file' is set to an open writable
        file-like object, then output is printed to that object instead.
        If 'file' is set to 'None', then no output is printed.
    'exclude_mutations' is an optional argument. If it has a value other
        than 'None', then it should be a dictionary. It is keyed by
        residue numbers (0, 1, 2, 3, ...). Values are lists of amino
        acid codes. Sequences that mutate indicated residues to amino
        acids in the lists for that residue are not allowed to be taken
        as steps on the walk.
    The returned variable is the 3-tuple '(seq, f, n)'.  'seq' is the 
	fitness of the best protein reached by the walk.  'f' is the
	fitness of 'seq'.  'n' is the number of steps used in the
	adaptive walk."""
    # Do some error checking on the input variables
    if not (fitness_cutoff == None or isinstance(fitness_cutoff, (int, float))):
        raise EvolutionError, "Invalid 'fitness_cutoff' of %r." % fitness_cutoff
    if not steepest_ascent in [True, False]:
        raise EvolutionError, "Invalid 'steepest_ascent' of %r." % steepest_ascent
    if not (isinstance(maxsteps, int) and maxsteps > 0):
        raise EvolutionError, "Invalid 'maxsteps' of %r." % maxsteps
    # make sure we can evaluate the fitness of the sequence
    try:
        initialf = fitness.Fitness(initialseq)
    except:
        raise EvolutionError, "Cannot evaluate the fitness of 'initialseq' with 'fitness'."
    # print initial information about walk
    if file != None:
        if steepest_ascent:
            file.write("Performing a steepest ascent walk of up to %d steps.\n" % maxsteps)
        else:
    	    file.write("Performing an adaptive walk of up to %d steps.\n" % maxsteps)
        fitness.Info(file)
        file.write("The initial sequence is %s, and has a fitness of %.3f.\n" % (''.join(initialseq), initialf))
        file.write("STEP\tFITNESS\tSEQUENCE\n")
        file.write("0\t%.5f\t%s\n" % (initialf, ''.join(initialseq)))
        file.flush()
    # begin the adaptive walk
    seq = list(initialseq)
    oldfit = initialf
    n = 0
    for istep in range(1, maxsteps + 1):
        if fitness_cutoff != None and oldfit >= fitness_cutoff:
            break # we have reached the fitness cutoff
        mutants = sequences.NMutants(seq, 1, 'ALL')
        if exclude_mutations:
            allowed = []
            for mutseq in mutants:
                for (r, aalist) in exclude_mutations.iteritems():
                    if mutseq[r] in aalist:
                        break
                else:
                    allowed.append(mutseq)
            mutants = allowed
        if not steepest_ascent: # shuffle the order of these mutants
            random.shuffle(mutants)
        bestseq = newfit = None
        for mutseq in mutants:
            f = fitness.Fitness(mutseq)
            if bestseq == None or f > newfit:
                newfit = f
                bestseq = mutseq
            if not steepest_ascent and newfit > oldfit:
                break 
        if newfit <= oldfit and (not newfit == oldfit == 0.0): # no improvement, walk ends
            break
        else: # go to the next step of the walk
            if newfit == oldfit == 0.0:
                seq = random.choice(mutants)
            else:
                seq = bestseq
                oldfit = newfit
        n = istep
        if file != None:
            file.write("%d\t%.5f\t%s\n" % (n, oldfit, ''.join(seq)))
    # return the variables
    return (''.join(seq), oldfit, n) 
#---------------------------------------------------------------------------
def GenerateRandomFoldingSequences(length, numtries, temp, walksteps, numcontacts, stability_cutoff = 0.0, steepest_ascent = False):
    """Generates a list of folding sequences.

    Call is: 'seqlist = GenerateRandomFoldingSequences(length, numseqs, 
	temp, walksteps, numcontacts, [stability_cutoff = 0.0,
	steepest_ascent = False])'
    This method generates random sequences that fold to a stable structure. 
	It does this by choosing a random starting sequence and performing
	an adaptive walk from that sequence.
    'length' is an integer > 0 specifying the length of the sequences generated.
    'numtries' is an integer > 0 specifying how many folding sequences we
	try to generate.  The number actually generated may be lower if some
	of the attempted adaptive walks fail to find a folding sequence.
    'temp' is the temperature at which the sequences are folded.
    'walksteps' is the maximum number of steps attempted in the adaptive walks.
    'numcontacts' is the number of contacts which the folded sequence must
	have.  If 'numcontacts' is 'None', then any conformation with any
	number of contacts is OK.  If 'numcontacts' is an integer > 0, then
	we require the folded sequence to have this number of contacts.
	This is done by randomly choosing one of the unique conformations
	with this many contacts, and then selecting for the ability of
	the sequences to fold to this conformation.
    'stability_cutoff' specifies the stability that the folding sequences
	must achieve.  They must have a free energy of folding <= 
	'stability_cutoff'.  The adaptive walks terminate when this
	cutoff is reached.  The default value is 0.0.
    'steepest_ascent' specifies whether the adaptive walk used to find the
	folding sequences is an adaptive walk or a steepest ascent walk.
	The default value is 'False', meaning it is an adaptive walk.  Set
	this argument to 'True' to get a steepest ascent walk.
    The returned variable is the list 'seqlist'.  This is simply a list of
	all of the folding sequences, as strings, that have been found."""
    # do some checking on input variables
    if not isinstance(stability_cutoff, (int, float)):
	raise EvolutionError, "Invalid 'stability_cutoff' of %r." % stability_cutoff
    if not (isinstance(length, int) and length > 0):
	raise EvolutionError, "Invalid 'length' of %r." % length
    if not (isinstance(numtries, int) and numtries > 0):
	raise EvolutionError, "Invalid 'numtries' of %r." % numtries
    if not (isinstance(walksteps, int) and walksteps > 0):
	raise EvolutionError, "Invalid 'walksteps' of %r." % walksteps
    if not (numcontacts == None or (isinstance(numcontacts, int) and numcontacts > 0)):
	raise EvolutionError, "Invalid 'numcontacts' of %r." % numcontacts 
    if not (isinstance(temp, (float, int)) and temp > 0):
	raise EvolutionError, "Invalid 'temp' of %r." % temp
    # create the 'Conformations' object to fold the proteins
    c = conformations.Conformations(length)
    if numcontacts == None:
	f = fitness.Fitness(temp, c, 'negstability', None)
    elif not c.UniqueConformations(numcontacts):
	raise EvolutionError, "No conformations with %r contacts." % numcontacts
    seqlist = []
    # begin trying to generate the conformations
    for itry in range(numtries):
	# recreate the fitness object if necessary
	if numcontacts:
	    f = fitness.Fitness(temp, c, 'negstability', random.choice(c.UniqueConformations(numcontacts)))
	# create the initial sequence
	s = sequences.RandomSequence(length)
	# do the adaptive walk
	(bestseq, bestfit, nsteps) = AdaptiveWalk(s, f, walksteps, fitness_cutoff = - stability_cutoff, steepest_ascent = steepest_ascent)
	# see if the sequence is stable enough
	if bestfit >= -stability_cutoff: 
	    seqlist.append(''.join(bestseq))
    # return the list of folding sequences
    return seqlist
#---------------------------------------------------------------------------
def FracViable(seq, fitness, fitness_cutoff, nmutations, nsampled = 'ALL'):
    """Computes the fraction of mutants that are viable.

    Call is: 'frac = FracViable(seq, fitness, fitness_cutoff, nmutations, 
	[nsampled = 'ALL'])'
    'seq' is the wild-type sequence which we are mutating, as a string or list.
    'fitness' is the object used to evaluate the fitness of the protein
	sequences.  It must have a method 'Fitness(seq)' such that
	'fitness.Fitness(seq)' returns a number representing the fitness
	of protein sequence 'seq'.
    'fitness_cutoff' specifies the minimum fitness a protein must have
	to be viable.  Any protein with fitness >= 'fitness_cutoff' is
	viable; otherwise it is not viable.
    'nmutations' specifies how many mutations to perform.  If it is one
	we sample single mutants, if it is two we sample double mutants,
	et cetera.  
    'nsampled' specifies how many of the mutants we sample.
	The default value is 'ALL', which means we sample all mutants.
	'ALL' is a valid value when 'nmutations' is one or two.  'nsampled'
	can also be an integer > 0, which means that we sample 'nsampled'
	randomly chosen mutants with 'nmutations' mutations.  In this case,
	the same mutation may be sampled multiple times.
    The returned variable is the number 'frac'.  This is simply the fraction
	of all of that mutants sampled that were viable."""
    # do some error checking on the initial variables
    if not isinstance(fitness_cutoff, (int, float)):
	raise EvolutionError, "Invalid 'fitness_cutoff' of %r." % fitness_cutoff
    # compute the wild-type fitness
    wt_fit = fitness.Fitness(seq)
    # compute the spectrum of mutation induced fitness changes
    spectrum = MutationSpectrum(seq, fitness, nmutations, nsampled)
    # compute the fraction viable
    nviable = 0
    for imut in spectrum:
	if wt_fit + imut >= fitness_cutoff:
	    nviable += 1
    # return the fraction
    try:
	return float(nviable) / len(spectrum)
    except ZeroDivisionError:
	raise EvolutionError, "No sequences for which to evaluate viability."
#----------------------------------------------------------------------
def MutationSpectrum(seq, fitness, nmutations, nsampled = 'ALL'):
    """Computes the spectrum of effects of mutations.

    Call is: 'spectrum = MutationSpectrum(seq, fitness, nmutations, 
	[nsampled = 'ALL'])'
    Creates a list of the effects of mutations on the fitness of a sequence.
    'seq' is the wild-type sequence which we are mutating, as a string or list.
    'fitness' is the object used to evaluate the fitness of the protein
	sequences.  It must have a method 'Fitness(seq)' such that
	'fitness.Fitness(seq)' returns a number representing the fitness
	of protein sequence 'seq'.
    'nmutations' specifies how many mutations to perform.  If it is one
	we sample single mutants, if it is two we sample double mutants,
	et cetera.  
    'nsampled' specifies how many of the mutants we sample.
	The default value is 'ALL', which means we sample all mutants.
	'ALL' is a valid value when 'nmutations' is one or two.  'nsampled'
	can also be an integer > 0, which means that we sample 'nsampled'
	randomly chosen mutants with 'nmutations' mutations.  In this case,
	the same mutation may be sampled multiple times.
    The returned variable is the list 'spectrum'.  Each entry in this
	list is the change in fitness, computed as the fitness of
	the mutant protein minus the fitness of the wild-type protein,
	for a particular mutation."""
    # do some error checking on the initial variables
    if not (isinstance(nmutations, int) and nmutations > 0):
	raise EvolutionError, "Invalid 'nmutations' of %r." % nmutations
    if nsampled == 'ALL' and nmutations in [1, 2]:
	pass
    elif not (isinstance(nsampled, int) and nsampled > 0):
	raise MutantError, "Invalid 'nsampled' of %r." % nsampled
    # make sure we can evaluate the fitness of the sequence
    try:
	wt_fit = fitness.Fitness(seq)
    except:
	error_msg = sys.exc_info()[1]
	raise EvolutionError, "Cannot evaluate the fitness of 'seq' with 'fitness': %r" % error_msg
    # now start sampling the mutants
    spectrum = []
    mutants = sequences.NMutants(seq, nmutations, nsampled)
    for mutseq in mutants:
	mut_fit = fitness.Fitness(mutseq)
	spectrum.append(mut_fit - wt_fit)
    return spectrum
#---------------------------------------------------------------------------
# End evolution.py
