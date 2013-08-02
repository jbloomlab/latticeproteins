#!/usr/bin/python
# Begin sequences.py
#---------------------------------------------------------------------------
"""Module for lattice protein sequences.

Written by Jesse Bloom, 2004."""
#---------------------------------------------------------------------------
import random, shelve, os
#---------------------------------------------------------------------------
class SequenceError(Exception):
    """Error with a lattice protein sequence."""
#---------------------------------------------------------------------------
# codes for all residues
_residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
assert len(_residues) == 20
#---------------------------------------------------------------------------
def MostAbundant(population):
    """Returns the most abundant sequence in a population.

    Call is: 'seq = MostAbundant(population)'
    'population' is a list of protein sequences.
    'seq' is returned as the sequence that is most abundant in 'population'.
	If there are several equally abundant sequences, just one of
	them is returned."""
    seq = population[0]
    n = population.count(seq)
    for seq2 in population[1 : ]:
	if seq2 != seq:
	    if population.count(seq2) > n:
		n = population.count(seq2)
		seq = seq2
    return seq
#---------------------------------------------------------------------------
def PairwiseHammingDistances(seqlist):
    """Computes the pairwise Hamming distances between many sequences.

    Call is: 'dlist = PairwiseHammingDistances(seqlist)'
    'seqlist' is a list of sequences all of the same length.
    'dlist' is returned as a list of numbers representing the
        Hamming distances between all pairs of sequences."""
    if not (isinstance(seqlist, list) and len(seqlist) > 1):
        raise SequenceError, "'seqlist' is not a list of at least 2 entries."
    length = len(seqlist[0])
    dlist = []
    for i1 in range(len(seqlist)):
        seq1 = seqlist[i1]
        if len(seq1) != length:
            raise SequenceError, "Invalid length sequence of %r." % seq1
        for i2 in range(i1 + 1, len(seqlist)):
            seq2 = seqlist[i2]
            if len(seq2) != length:
                raise SequenceError, "Invalid length sequence of %r." % seq1
            dlist.append(HammingDistance(seq1, seq2))
    if len(dlist) != len(seqlist) * (len(seqlist) - 1) / 2:
        raise SequenceError, "Incorrect number of distances."
    return dlist
#---------------------------------------------------------------------------
def HammingDistance(seq1, seq2):
    """Returns the Hamming distance between two sequences.

    Call is: 'd = HammingDistance(seq1, seq2)'
    'seq1' and 'seq2' are two sequences of the same length.
    'd' is returned as the Hamming distance between these two sequences."""
    if len(seq1) != len(seq2):
	raise SequenceError, "Sequences differ in length."
    d = 0
    for i in range(len(seq1)):
	if seq1[i] != seq2[i]:
	    d += 1
    return d
#---------------------------------------------------------------------------
def RandomSequence(length):
    """Returns a random sequence of the specified length.

    Call is 's = RandomSequence(length)'
    'length' is an integer >= 1.  Returns a sequence of length 'length'
	as a list of randomly chosen residues."""
    if not (isinstance(length, int) and length > 0):
	raise SequenceError, "Invalid sequence length of %r." % length
    s = [random.choice(_residues) for i in range(length)]
    return s
#------------------------------------------------------------------------
def MutateSequence(seq, mutrate):
    """Mutates a protein sequence.

    Call is: 'seqnew = MutateSequence(seq, mutrate)'
    'seq' is a protein sequence, specified as either a string or a list.
    Mutates each residue in 'seq' to some different residue with
	probability 'mutrate'.  So 'mutrate' is the per residue 
        mutation rate.
    Returns the new sequence as the list 'seqnew'."""
    mutated = False
    for ires in range(len(seq)):
	if random.random() < mutrate:
            if not mutated:
                mutated = True
                newseq = list(seq)
	    newres = random.choice(_residues)
	    while newres == seq[ires]:
		newres = random.choice(_residues)
            newseq[ires] = newres
    if mutated:
        return newseq
    else:
        return seq
#------------------------------------------------------------------------
def NMutants(seq, nmutations, nsequences):
    """Returns sequences with a specified number of mutations.

    Call is: 'seqlist = NMutants(seq, nmutations, nsequences)'
    'seq' is a string or list specifying the protein we wish to mutate.
    'nmutations' is the number of mutations each mutant of 'seq' should
	have.  It must be <= 'len(seq)' and > 0.
    'nsequences' is the number of mutant sequences to make.  It can be
	'ALL', in which case we make all possible mutants with 'nmutations',
	or it can be some positive integer in which case we make this
	many randomly chosen mutants with 'nmutations' mutations.
	'ALL' is only a valid option only when 'nmutations' is 1 or 2.
    The mutant sequences are returned in the list 'seqlist'.  Each entry
	in 'seqlist' is a list representing a mutant sequence."""
    if not (0 < nmutations <= len(seq)):
	raise SequenceError, "Invalid 'nmutations' of %r." % nmutations
    seqlist = []
    if nsequences == 'ALL':
	if nmutations == 1:
	    for ires in range(len(seq)):
		for mutres in _residues:
		    if mutres != seq[ires]:
			newseq = list(seq)
			newseq[ires] = mutres
			seqlist.append(newseq)
	elif nmutations == 2:
	    for ires in range(len(seq)):
		for imutres in _residues:
		    if imutres != seq[ires]:
			for jres in range(ires + 1, len(seq)):
			    for jmutres in _residues:
				if jmutres != seq[jres]:
				    newseq = list(seq)
				    newseq[ires] = imutres
				    newseq[jres] = jmutres
				    seqlist.append(newseq)
	else:
	    raise SequenceError, "'nsequences' cannot be 'ALL' when 'nmutations' is %r." % nmutations
    elif isinstance(nsequences, int) and nsequences > 0:
	for imutant in range(nsequences):
	    newseq = list(seq)
	    for imut in range(nmutations):
		ires = random.choice(range(len(seq)))
		while newseq[ires] != seq[ires]:
		    ires = random.choice(range(len(seq)))
	      	mutres = random.choice(_residues)
    		while mutres == seq[ires]:
       		    mutres = random.choice(_residues)
    		newseq[ires] = mutres
	    seqlist.append(newseq)
    else:
	raise SequenceError, "Invalid 'nsequences' of %r." % nsequences
    return seqlist
#---------------------------------------------------------------------------
# '_foldingsequences_database' is a database that stores folding sequences.
# It is accessed by 'SaveFoldingSequence' and 'GetFoldingSequence'.
# The database is keyed by string conversions of integer sequence lengths.
# The values of the sequence length keys are dictionaries keyed by temperatures.
# The values of the temperature keys are dictionaries keyed by number of contacts.
# The values of the contact keys are lists of 3-tuples giving the free energy of 
# folding, the sequence, and the conformation for all sequences.
_foldingsequences_database = '_foldingsequences.database'
#---------------------------------------------------------------------------
def SaveFoldingSequence(seq, temp, dGf, numcontacts, conf):
    """Saves a folding sequence in the database of folding sequences.

    Call is: 'SaveFolding(seq, temp, dGf, numcontacts, conf)'
    The sequence is saved in the database of folding sequences, where it
	can be accessed with 'GetFoldingSequence()'.
    'seq' is a protein sequence that folds to conformation 'conf'
	at temperature 'temp' with a free energy of folding of 'dGf'.
	The number of contacts for conformation 'conf' is 'numcontacts'."""
    # Do some error checking on the input variables
    if len(seq) != len(conf) + 1:
	raise SequenceError, "'seq' and 'conf' have incompatible lengths."
    if not (isinstance(numcontacts, int) and numcontacts > 0):
	raise SequenceError, "Invalid 'numcontacts' of %r." % numcontacts
    if not (isinstance(temp, float) and temp > 0):
	raise SequenceError, "Invalid 'temp' of %r." % temp
    # save in the database
    try:
	# open the databae
	database = shelve.open(_foldingsequences_database, 'c')
	# access the relevant sequence list
	dbkey = str(len(seq)) # key into the database, sequence length as string
	try:
	    length_dict = database[dbkey]
	except KeyError:
	    length_dict = {}
	try:
	    temp_dict = length_dict[temp]
	except KeyError:
	    temp_dict = {}
	    length_dict[temp] = temp_dict
	try:
	    seq_list = temp_dict[numcontacts]
	except KeyError:
	    seq_list = []
	    temp_dict[numcontacts] = seq_list
	# add the entry to the sequence list
	seq_list.append((dGf, seq, conf))
	# store the updated information
	database[dbkey] = length_dict
    finally:
	database.close()
#---------------------------------------------------------------------------
def GetFoldingSequence(length, temp, dGf, numcontacts = None, conf = None, nget = 1):
    """Gets a folding sequence from the database of folding sequences.

    Call is: 's = GetFoldingSequence(length, temp, dGf, 
	[numcontacts = None, conf = None, nget = 1])'
    This function tries to find a protein in the database of folding sequence 
	(entries added by 'SaveFoldingSequence') meeting the specified 
	criteria, and returns its sequence.  If no such protein can be 
	found, it returns 'None'.
    'length' is an integer specifying the lengths of sequences we consider.
    'temp' can be a number, in which case we consider sequences that
	meet the folding criterium specified by 'dGf' at temperature 
	'temp'.  Or it can be a 2-tuple, '(mint, maxt)', in which 
	case we consider sequences that meet the folding criterium
	specified by 'dGf' at some temperature t: mint <= t <= maxt.
    'dGf' can be a number, in which case we consider any sequences
	with folding free energies <= 'dGf'.  Or it can be a 2-tuple,
	'(ming, maxg)', in which case we consider sequences with 
	folding free energies g: ming <= g <= maxg.
    'numcontacts' is on optional argument with a default value of 'None'.  If
	it has the value of 'None', then no selection based on the number
	of contacts is performed.  If it is set to another value, it
	can be an integer, in which case we consider only sequences
	that fold to structures with 'numcontacts' contacts.  Or it can
	be a 2-tuple of integers, '(minn, maxn)', in which case we consider
	only sequences that fold to structures that fold with a number
	of contacts n: minn <= n <= maxn.
    'conf' is an optional argument with a default value of 'None'.  If it has
	the value of 'None', then no selection based on conformation is
	performed.  If it is set to another value, it should be set to
	a string indicating a conformation for a protein of length
	'length'.  In this case, only sequences with conformations of
	'conf' are considered as candidates for return.
    'nget' is an optional argument specifying how many sequences to return.
	If it has its default value of 1, a single sequence is returned, which
	is chosen at random from all of the sequences in the database
	that meet the specified criteria.  If 'nget' is an integer > 1,
	a list of sequences is returned.  This list contains 'nget'
	sequences that meet the specified criteria, or if there are less
	than 'nget' such sequences, returns all that are present.  If 'nget'
	is 'ALL', returns all sequences in the database that meet the 
	specified criteria.
    's' is returned as 'None' if no sequences can be found that meet the
	criteria.  If 'nget' is 1, it is a sequence returned as either a
	list or string.  If 'nget' is greater than 1, it is a list
	of sequences (even if there is only one sequence in the list)."""
    # if there is no database, just return 'None'
    if not os.path.isfile(_foldingsequences_database):
	return
    candidates = [] # a list of all sequences meeting the criteria
    # get the dictionary for this sequence list from the dictionary
    try:
	database = shelve.open(_foldingsequences_database)
	length_dict = database[str(length)]
    finally:
	database.close()
    # now search for sequences in the this dictionary
    # look for candidates with the right temperature
    temp_dict_list = []
    if isinstance(temp, (int, float)):
	try:
    	    temp_dict_list.append(length_dict[temp])
	except KeyError:
	    pass
    elif isinstance(temp, tuple) and len(temp) == 2:
	for (t, d) in length_dict.iteritems():
	    if temp[0] <= t <= temp[1]:
		temp_dict_list.append(t)
    else:
	raise SequenceError, "Invalid 'temp' of %r." % temp
    # look for candidates with the right number of contacts
    numcontacts_list = []
    for d in temp_dict_list:
	if numcontacts == None:
	    for (n, l) in d.iteritems():
		numcontacts_list += l # consider all numbers of contacts
	elif isinstance(numcontacts, int):
	    try:
    		numcontacts_list += d[numcontacts]
	    except KeyError:
		pass
	elif isinstance(numcontacts, tuple) and len(numcontacts) == 2:
	    for (n, l) in d.iteritems():
		if numcontacts[0] <= n <= numcontacts[1]:
		    numcontacts_list += l
	else:
	    raise SequenceError, "Invalid 'numcontacts' of %r." % numcontacts
    # look for sequences with the correct dGf
    if isinstance(dGf, (int, float)):
    	for (g, seq, c) in numcontacts_list:
	    if g <= dGf:
		if conf == None:
    		    candidates.append(seq)
		elif isinstance(conf, str) and len(conf) == length - 1:
		    if conf == c:
			candidates.append(seq)
		else:
		    raise SequenceError, "Invalid 'conf' of %r." % conf
    elif isinstance(dGf, tuple) and len(dGf) == 2:
    	for (g, seq, c) in numcontacts_list:
	    if dGf[0] <= g <= dGf[1]:
		if conf == None:
    		    candidates.append(seq)
		elif isinstance(conf, str) and len(conf) == length - 1:
		    if conf == c:
			candidates.append(seq)
		else:
		    raise SequenceError, "Invalid 'conf' of %r." % conf
    else:
    	raise SequenceError, "Invalid 'dGf' of %r." % dGf
    # if there are no sequences, return 'None'
    if not candidates:
	return
    # return the correct number of sequences
    if nget == 'ALL':
	return candidates
    elif nget == 1:
	return random.choice(candidates)
    elif isinstance(nget, int) and nget > 1:
	if nget >= len(candidates):
	    return candidates
	else:
	    random.shuffle(candidates)
	    return candidates[0 : nget]
    else:
	raise SequenceError, "Invalid 'nget' of %r." % nget
#---------------------------------------------------------------------------
# End sequences.py
