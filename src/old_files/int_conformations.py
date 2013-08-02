#!/usr/bin/python
# Begin conformations.py
#---------------------------------------------------------------------------
"""Module for finding and storing the conformations of a 2D lattice protein.

Written by Jesse Bloom, 2004."""
#-----------------------------------------------------------------------
import math, shelve, interactions, sys
#----------------------------------------------------------------------
class ConformationsError(Exception):
    """Error finding or storing a conformation."""
#----------------------------------------------------------------------
# '_loop_in_C' is a Boolean switch specifying whether we try to speed
# up the energy calculations by doing the looping with the C-extension
# 'contactlooper'.  'True' means we try to do this.
_loop_in_C = True 
if _loop_in_C:
    import contactlooper
#----------------------------------------------------------------------
class Conformations(object):
    """A class for storing all of the conformations of a lattice protein."""
    #-------------------------------------------------------------------
    # '_existing_instance' is a flag indicating whether we are returning an
    # existing instance or creating a new one
    _existing_instance = False 
    #------------------------------------------------------------------
    def __new__(cls, length, database_file):
        """Gets existing or creates new listing of conformations.

        Arguments are as to '__init__'.  If the specified 'Conformations'
            object is already present in the class's database returns it, 
            otherwise generates a new one."""
        try:
            database = shelve.open(database_file, 'c', protocol=2)
            if database.has_key(str(length)):
                c = database[str(length)]
                cls._existing_instance = True
                if not isinstance(c, cls):
                    raise ConformationsError, "'Conformation' read from database is invalid."
            else:
                cls._existing_instance = False 
                c = object.__new__(cls)
        finally:
            database.close()
        return c
    #------------------------------------------------------------------
    def __init__(self, length, database_file):
        """Creates a list of conformations for a protein of specified length.

        Call is: 'c = Conformations(length, database_file)
        'length' is an integer specifying the length of the protein
            for which we are computing the contacts.  It must be >= 2.
        'database_file' specifies the name of the database storing
            existing conformations.  If the conformation instance
            already exists in this database we return the
            existing instance, and if it doesn't we store it in the database.
        The created 'Conformations' object 'c' stores the contact
            lists and the number of conformations with these contact sets
            for all self-avoiding walks of length 'length'.  It can then 
            be used to compute the free energy of a protein folding to
            the lowest energy conformation."""
        # If '__new__' has returned an existing instance, just return
        if self._existing_instance:
            # sort the contact set information so that the conformations
            # with the most contacts come first
            n = len(self._contactsets)
            decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(n)]
            decorated_list.sort()
            decorated_list.reverse()
            self._contactsets = [decorated_list[i][1] for i in range(n)]
            self._contactsetdegeneracy = [decorated_list[i][2] for i in range(n)]
            self._contactsetconformation = [decorated_list[i][3] for i in range(n)]
            return
        # A listing of all class variables
        # 'self._length' is the length of the protein (self-avoiding walk)
        self._length = length 
        # there are 'self._numconformations[i]' conformations with 'i' contacts
        self._numconformations = {}
        # 'self._contactsets' is a list of contact sets.  'self._contactsets[i]'
        # is the contact set for contact i.  It is a list of numbers.
        # 'x = self._contactsets[i]' describes the residues in contact
        # in contact 'i'.  If this contact is between residues 'ires'
        # and 'jres', then 'x = self._length * ires + jres' where
        # 0 <= ires, jres < 'self._length', and ires < jres + 1
        self._contactsets = []
        # 'self._contactsetdegeneracy' is a list of integers giving the
        # degeneracy of the contact sets (the number of different conformations
        # with this contact set).  'self_contactsetdegeneracy[i]' is the
        # degeneracy of the contact set 'self._contactsets[i]'
        self._contactsetdegeneracy = []
        # 'self._contactsetconformation' is a list of the conformations
        # associated with each contact set.  If contact set 'self._contactsets[i]'
        # is degenerate ('self._contactsetdegeneracy[i]' > 1), the value
        # 'self._contactsetconformation[i]' is 'None'.  Otherwise, it
        # is the string representing the conformation that gives rise
        # to contact set 'self._contactsets[i]'.  The conformations are given 
        # such that 'self._contactsetconformation[i][j]'
        # gives the conformation of bond 'j' (0 <= j < 'self._length' - 1)
        # as 'U' (Up), 'R' (Right), 'D' (Down), or 'L' (Left). 
        # We require the first bond to be Up, and the first
        # non-Up bond to be Right. 
        self._contactsetconformation = []
        # 'self._numcontactsets[i]' holds the number of different contact 
        # sets with 'i' contacts.
        self._numcontactsets = {}
        # 'self._foldedsequences = {}' holds recently folded sequences.  
        # The keys are the arguments to 'FoldSequence', and the items
        # are how many times the saved sequence was accessed and the 
        # information about the folded sequence.  Sequences in the keys
        # are strings.
        self._foldedsequences = {}
        #
        # Now do some error checking on input variables
        if not (isinstance(self._length, int) and self._length >= 2):
            raise ConformationsError, "Invalid 'length' of %r." % self._length
        #
        # The initial generation of conformations uses the variable
        # 'contactsets' which is a dictionary as described below.  Once
        # we have constructed this dictionary, we use it to create
        # 'self._contactsets', 'self._contactsetdegeneracy', and
        # 'self._contactsetconformation'.  The format of 'contactsets':
        # The keys are string representing the contact sets.  These
        # strings are the numbers that will appear in 'self._contactsets',
        # separated by spaces.
        # The values of 'contactsets' provides information
        # about the conformations encoding the contact sets.  If there is
        # a single conformation encoding a contact set 'cs', then
        # 'contactsets[cs]' is a string.  If the contact set stored as
        # 'contactsets[cs]' is coded for by multiple conformations,
        # then 'self._contactsets[cs]' is an integer > 1 that represents
        # the number of conformations coding for this contact set.
        contactsets = {}
        # 
        # Now being looping over all possible conformations
        # The initial conformation is all 'U'
        dx = {'U':0, 'R':1, 'D':0, 'L':-1}
        dy = {'U':1, 'R':0, 'D':-1, 'L':0}
        next = {'U':'R', 'R':'D', 'D':'L', 'L':'U'} 
        opposite = {'U':'D', 'D':'U', 'R':'L', 'L':'R'}
        n = self._length - 2 # index of last bond in 'conformation' 
        conformation = ['U' for i in range(n + 1)] 
        first_R = n # index of the first 'R' in the conformation
        ncount = 0
        while True: # keep finding new conformations
            # See if the current conformation has overlap
            # If no overlap, store the contact set
            x = y = j = 0
            res_positions = {(x, y):j} # keyed by coords, items are residue numbers
            res_coords = [(x, y)] # 'res_coords[j]' is coords of residue 'j'
            for c in conformation:
                x += dx[c]
                y += dy[c]
                coord = (x, y)
                if coord in res_positions: # overlap
                    # increment at the step that gave the problem
                    for k in range(j + 1, n + 1):
                        conformation[k] = 'U'
                    conformation[j] = next[conformation[j]]
                    while conformation[j] == 'U':
                        j -= 1
                        conformation[j] = next[conformation[j]]
                    if j == first_R and conformation[j] not in ['R', 'U']:
                        first_R -= 1
                        conformation[first_R] = 'R'
                        for k in range(j, n + 1):
                            conformation[k] = 'U'
                    break
                j += 1
                res_positions[coord] = j
                res_coords.append(coord)
            else: # loop finishes normally, this is a valid conformation
                # generate the contact set
                cs = []
                numcontacts = 0
                for j in range(len(res_coords)):
                    (x, y) = res_coords[j]
                    partners_list = []
                    for c in ['U', 'R', 'D', 'L']:
                        try:
                            k = res_positions[(x + dx[c], y + dy[c])]
                            if k > j + 1:
                                partners_list.append(k)
                        except KeyError:
                            pass
                    partners_list.sort()
                    jtimeslength = j * self._length
                    for k in partners_list:
                        cs.append("%d " % (jtimeslength + k))
                        numcontacts += 1
                cs = ''.join(cs)
                # add conformation to count
                try:
                    self._numconformations[numcontacts] += 1
                except KeyError:
                    self._numconformations[numcontacts] = 1
                # store contact set
                try:
                    contactsets[cs] += 1
                except KeyError:
                    contactsets[cs] = ''.join(conformation)
                    try:
                        self._numcontactsets[numcontacts] += 1
                    except KeyError:
                        self._numcontactsets[numcontacts] = 1 
                except TypeError:
                    contactsets[cs] = 2
                # generate the next conformation
                i = n
                conformation[i] = next[conformation[i]]
                while conformation[i] == 'U':
                    i -= 1
                    conformation[i] = next[conformation[i]]
                # make sure first non-'U' is 'R'
                if i == first_R and conformation[i] not in ['R', 'U']:
                    first_R -= 1
                    conformation[first_R] = 'R'
                    for j in range(i, n + 1):
                        conformation[j] = 'U'
            # see if we are done
            if first_R == 0:
                break
        #
        # Now use 'contactsets' to generate 'self._contactsets', 
        # 'self._contactsetdegeneracy', and 'self._contactsetconformation'
        for (cs, n_or_conf) in contactsets.iteritems():
            # convert 'cs' to the appropriate list
            clist = [int(x) for x in cs.split()]
            self._contactsets.append(clist)
            if isinstance(n_or_conf, str):
                self._contactsetdegeneracy.append(1)
                self._contactsetconformation.append(n_or_conf)
            else:
                self._contactsetdegeneracy.append(n_or_conf)
                self._contactsetconformation.append(None)
        #
        del contactsets
        # 
        # to make the energy evaluations quicker, sort so that
        # contact sets with more conformations are first:
        decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(len(self._contactsets))]
        del self._contactsets, self._contactsetdegeneracy, self._contactsetconformation
        decorated_list.sort()
        decorated_list.reverse()
        self._contactsets = [decorated_list[i][1] for i in range(len(decorated_list))]
        self._contactsetdegeneracy = [decorated_list[i][2] for i in range(len(decorated_list))]
        self._contactsetconformation = [decorated_list[i][3] for i in range(len(decorated_list))]
        #
        # store the contact set in the database
        try:
            database = shelve.open(database_file, 'c', protocol=2)
            database[str(self._length)] = self
        finally:
            database.close()
    #------------------------------------------------------------------
    def Length(self):
        """Returns the length of the protein these conformations are for.

        Call is: 'length = c.Length()'."""
        return self._length
    #------------------------------------------------------------------
    def FoldSequence(self, seq, temp, target_conf = None, numsaved = 1e6, dGf_cutoff = None):
        """Folds a protein sequence.

        Call is: '(dGf, conf, numcontacts) = c.FoldSequence(seq, temp,  
            [target_conf = None, numsaved = 1e6])'
        'seq' is the sequence of the protein to be folded as one-letter amino
            acid codes.  It should be a string or list of length 'c.Length()'.
        'temp' is the temperature at which the protein is to be folded.  It
            must be a number > 0.  It represents a reduced temperature, scaled
            so that a value of 1 represents 273 K.
        'target_conf' is an optional argument specifying the target
            conformation to which we fold the protein.  If it has its
            default value of 'None', we fold the protein to its lowest
            energy conformation.  If 'target_conf' is not 'None', it
            must be set to a valid conformation string (as described below
            in the description of 'conf') describing a unique conformation
            (only one conformation with this contact set.  The returned free 
            energy 'dGf' is then for folding to 'target_conf', and 'conf'
            is the same as 'target_conf'.
        'dGf_cutoff' is meaningful iff 'target_conf' is not 'None'.
            In that case, we terminate computing dGf once we can
            ensure that it is greater than dGf_cutoff.  In this
            case, dGf is returned as 'None'.  This option only
            works when we are looping in C (hardcoded constant at
            the beginning of this module).
        'numsaved' is an option utilizing the fact that the method can
            save the information for recently folded sequences instead
            of refolding them.  At any time, up to the last 'numsaved'
            unique sequences that were folded under different conditions
            are saved.  If 'seq' folded with the options set by 
            'temp' and 'target_conf'is one of these
            saved sequences, it can be taken from the saved sequences
            and is not refolded.  As soon as more than 'numsaved'
            sequences have been folded since the last clearing of
            the folded sequences, half of the saved sequences are
            deleted from the list of saved sequences.  The sequences
            that are deleted are the ones that have been accessed
            the least frequenty since the last deletion of saved sequences.
        The returned 3-tuple specifies the free energy of folding to the
            lowest energy conformation, the conformation if there
            is a single unique lowest energy conformation, and the
            number of contacts in that conformation.  'dGf'
            is the free energy of folding at temperature 'temp', computed
            from the partition function.  If there is a single
            lowest energy conformation, then 'conf' is a string
            specifying that conformation as follows: the string consists
            of the characters 'U', 'D', 'L', or 'R'.  'conf' is of length
            'len(seq)' - 1, and 'conf[j]' is the conformation of bond 'j'
            (the bond between residues j and j + 1).  'U' means Up, 
            'D' means Down, 'R' means Right, and 'L' means Left.
            'numcontacts' is the number of contacts in the lowest
            energy conformation.  If there is no single lowest energy
            conformation, the 'conf' is 'None'."""
        # if we have more saved sequences than 'numsaved', get rid of them
        if numsaved < len(self._foldedsequences):
            decorateditems = []
            for i in self._foldedsequences.iteritems():
                decorateditems.append((i[1][0], i))
            self._foldedsequences = {}
            decorateditems.sort()
            decorateditems.reverse()
            for i in decorateditems[0 : int(numsaved / 2)]:
                self._foldedsequences[i[1][0]] = (0, i[1][1][1])
        # Get the sequence if it stored
        savekey = (''.join(seq), temp, target_conf, dGf_cutoff)
        try:
            self._foldedsequences[savekey] = (self._foldedsequences[savekey][0] + 1, self._foldedsequences[savekey][1])
            return self._foldedsequences[savekey][1]
        except KeyError:
            pass
        # do some error checking on the input variables
        if len(seq) != self._length:
            raise ConformationsError, "Invalid 'seq' length: is %r but should be %r." % (len(seq), self._length)
        if target_conf != None:
            if not (isinstance(target_conf, str) and len(target_conf) == len(seq) - 1):
                raise ConformationsError, "Invalid 'target_conf' of %r." % target_conf
                if dGf_cutoff != None:
                    if not isinstance(dGf_cutoff, (int, float)):
                        raise ConformationsError, "Invalid 'dGf_cutoff' of %s." % dGf_cutoff
        try:
            temp = float(temp)
            if temp <= 0.0:
                raise ConformationsError, "'temp' is <= 0: %r." % temp
        except KeyError:
            raise ConformationsError, "Invalid 'temp' of %r." % temp
        # create 'res_interactions'.  'res_interactions[j]' holds the energy
        # for the interaction 'j' as specified in 'self._contactsets[i][j]'
        res_interactions = [0.0 for i in range(self._length**2)]
        for ires in range(self._length):
            itimeslength = ires * self._length
            for jres in range(ires + 1, self._length):
                respair = "%s%s" % (seq[ires], seq[jres])
                res_interactions[itimeslength + jres] = interactions.miyazawa_jernigan[respair]  
        # now loop over all contact sets
        partitionsum = 0.0 # the partition sum
        contactsets = self._contactsets
        contactsetdegeneracy = self._contactsetdegeneracy
        contactsetconformation = self._contactsetconformation
        # we write out separate loops for the two possible value of 'target_conf'
        # first for the case where 'target_conf' is 'None':
        if target_conf == None:
            if _loop_in_C: # use the fast 'contactlooper' C-extension
                (minE, ibest, partitionsum) = contactlooper.NoTargetLooper(res_interactions, contactsets, contactsetdegeneracy, float(temp))
            else: # do the looping in python
                # initially set minimum to the first contact set:
                minE = 0.0
                for pair in contactsets[0]:
                    minE += res_interactions[pair]
                ibest = 0
                for i in range(len(self._contactsets)): 
                    e_contactset = 0.0 # energy of this contact set
                    # loop over all contact pairs in the contact set
                    for pair in contactsets[i]:
                        e_contactset += res_interactions[pair]
                        # add the energy for this contact set to the partition sum
                        partitionsum += math.exp(-e_contactset / temp) * contactsetdegeneracy[i]
                    if e_contactset < minE:
                        minE = e_contactset
                        ibest = i
            conf = contactsetconformation[ibest]
            numcontacts = len(contactsets[ibest])
        # now for the case where 'target_conf' is a conformation
        else:
            if _loop_in_C: # use the fast 'contactlooper' C-extension
                if dGf_cutoff != None:
                    (minE, ibest, partitionsum) = contactlooper.TargetLooper(res_interactions, contactsets, contactsetdegeneracy, float(temp), target_conf, contactsetconformation, 1, float(dGf_cutoff))
                else:
                    (minE, ibest, partitionsum) = contactlooper.TargetLooper(res_interactions, contactsets, contactsetdegeneracy, float(temp), target_conf, contactsetconformation, 0, 0.0)
                numcontacts = len(contactsets[ibest])
            else: # do the looping in python
                minE = conf = numcontacts = None # lowest energy sequence properties 
                for i in range(len(self._contactsets)): 
                    e_contactset = 0.0 # energy of this contact set
                    # loop over all contact pairs in the contact set
                    for pair in contactsets[i]:
                        e_contactset += res_interactions[pair]
                        # add the energy for this contact set to the partition sum
                        partitionsum += math.exp(-e_contactset / temp) * contactsetdegeneracy[i]
                    if target_conf == contactsetconformation[i]:
                        minE = e_contactset
                        numcontacts = len(contactsets[i]) 
            if minE == None:
                raise ConformationsError, "'target_conf' is not a unique conformation."
            conf = target_conf
        if dGf_cutoff != None and target_conf and partitionsum < 0:
            dGf = None
        else:
            gu = - temp * math.log(partitionsum - math.exp(-minE / temp)) 
            dGf = minE - gu
        # store
        returnkey = (dGf, conf, numcontacts)
        self._foldedsequences[savekey] = (0, returnkey)
        # return the variables
        return returnkey 
    #------------------------------------------------------------------
    def UniqueConformations(self, numcontacts):
        """Gets all unique conformations with specified number of contacts.

        Call is: 'clist = c.UniqueConformations(numcontacts)'
        'numcontacts' is an integer >= 0.  
        The returned list 'clist' is of all unique conformations
            with exactly 'numcontacts' contacts.  A conformation
            is "unique" if it is the only conformation that gives
            rise to its particular contact set.  If there are
            no unique conformations with 'numcontacts' contacts,
            'clist' is an empty list.  Conformations are specified
            as strings of 'U', 'R', 'L', and 'D' as described in 
            'FoldSequence'."""
        if not (isinstance(numcontacts, int) and numcontacts >= 0):
            raise ConformationsError, "Invalid 'numcontacts' of %r." % numcontacts
        clist = []
        for i in range(len(self._contactsets)):
            if self._contactsetdegeneracy[i] == 1:
                if len(self._contactsets[i]) == numcontacts:
                    clist.append(self._contactsetconformation[i])
        return clist
    #------------------------------------------------------------------
    def MaxContacts(self):
        """Gets the most contacts of any conformation.

        Call is: 'n = c.MaxContacts()'
        'n' is returned as the number of contacts for the conformation
            with the most contacts."""
        n = 0
        for cs in self._contactsets:
            if len(cs) > n:
                n = len(cs)
        return n 
    #------------------------------------------------------------------
    def NumConformations(self, contacts = None):
        """Returns the number of conformations.

        Call is: 'n = c.NumConformations([contacts = None])'
        If 'contacts' has its default value of 'None', returns the total
            number of conformations (self-avoiding walks).
        If 'contacts' has an integer value, returns the number of conformations
            with 'contacts' contacts.  If there are no walks with this number
            of contacts, returns 0."""
        if contacts == None:
            n = 0
            for x in self._numconformations.itervalues():
                n += x
        else:
            try:
                n = self._numconformations[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise ConformationsError, "Invalid 'contacts' of %r." % contacts
        return n
    #------------------------------------------------------------------
    def NumContactSets(self, contacts = None):
        """Returns the number of unique contact sets.

        Call is: 'n = c.ContactSets([contacts = None])'
        If 'contacts' has its default value of 'None', returns the total
            number of unique contact sets (defined as the list of all
            contacts of non-adjacent residues).
        If 'contacts' has an integer value, returns the number of unique
            contact sets with 'contacts' contacts.  If there are no
            contact sets with this number of contacts, returns 0."""
        if contacts == None:
            return len(self._contactsets)
        else:
            try:
                return self._numcontactsets[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise Conformationserror, "Invalid 'contacts' of %r." % contacts
#---------------------------------------------------------------------------
_saved_combinations = {} # saved combinations for 'BindLigand'.
# The keys are the 2-tuples '(protconf, ligandconf)' and the items
# are how many times the saved sequence was accessed and the 
# information about the folded sequence.  
_saved_exactcombinations = {} # saved exact combinations for 'BindLigand'.
# The keys are the 4-tuples '(prot, protconf, ligand, ligandconf)' and
# the items are how many times the combination was accessed and the
# information about the best binding position.
#----------------------------------------------------------------------
def BindLigand(prot, protconf, ligand, ligandconf, numsaved = 5000, numsavedexact = 5000):
    """Finds the lowest energy for a ligand binding to a lattice protein.

    Call is: '(be, xshift, yshift, conf) = BindLigand(prot, protconf, ligand, 
        ligandconf, [numsaved = 5000, numsavedexact = 5000])'
    'prot' is the string giving the sequence of the protein to which the ligand
        is being bound.
    'protconf' is a string specifying the conformation of 'prot'.  The string
        consists of the characters 'U', 'D', 'L', or 'R', and is of
        length 'len(protseq)' - 1.  'protconf[j]' is the conformation of bond
        'j' (the bond between residues j and j + 1).  'U' means Up,
        'D' means Down, 'R' means Right, and 'L' means Left.
    'ligand' is the string giving the sequence of the ligand which is 
        being bound to the protein.
    'ligandconf' is a string specifying the conformation of 'ligand', in
        the same format as 'protconf'.
    'be' is returned as the binding energy of 'ligand' to 'prot' in
        the lowest energy binding position.  This is computed by searching
        over all translational and rotational possible positions of the
        ligand relative to the protein.
    'xshift', 'yshift', and 'conf' specify the lowest energy binding position.
        To compute the position of the ligand relative to the protein,
        start at the position of the first residue of the protein.  Move
        'xshift' lattice positions to the right, and 'yshift' lattice positions
        up. Then construct the ligand using the ligand sequence given by 'ligand'
        and the conformation given by 'conf', with the first ligand residue
        in the position specified by 'xshift' and 'yshift'.  The reason that 
        'conf' might be different from 'ligandconf' is that the ligand might
        be rotated.
    'numsaved' specifies that we save the possible binding positions 
        for a given ligand conformation and protein conformation.  This
        means that we do not have to regenerate these positions repeatedly
        if the function is called many times with the same protein 
        conformation and ligand conformation.  As soon as more than 'numsaved'
        protein/ligand combinations have been analyzed since the last clearing 
        of the positions, half of the saved positions are deleted from the
        record of saved sequences.  The combinations that are deleted are
        the ones that have been accessed the least frequently since the
        last deletion of combinations.
    'numsavedexact' specifies that we save the binding energies for a 
        given protein and ligand sequence and conformation combination.
        This means that in addition to not having to regenerate the
        positions repeatedly (ensured by 'numsaved'), we also do not
        have to recalculate the best binding position.  As soon as more
        than 'numsavedexact' combinations have been analyzed since the last 
        clearing of the positions, half of the saved positions are deleted
        from the record of saved combinations.  The combinations that 
        are deleted are the ones that have been accessed the least
        frequently since the last deletion of combinations."""
    # First see if we have this exact combination saved
    # If we have more saved exact combinations than 'numsavedexact',
    # get rid of half of them.
    if numsavedexact < len(_saved_exactcombinations):
        decorateditems = [(i[1][0], i) for i in _saved_exactcombinations.iteritems()]
        decorateditems.sort()
        decorateditems.reverse()
        for i in decorateditems[int(numsaved / 2) : ]:
            del _saved_exactcombinations[i[1][0]]
        for i in decorateditems[0 : int(numsaved / 2)]:
            _saved_exactcombinations[i[1][0]] = (0, i[1][1][1])
    exactsavekey = (''.join(prot), ''.join(protconf), ''.join(ligand), ''.join(ligandconf))
    try:
        _saved_exactcombinations[exactsavekey] = (_saved_exactcombinations[exactsavekey][0] + 1, _saved_exactcombinations[exactsavekey][1])
        return _saved_exactcombinations[exactsavekey][1]
    except KeyError:
        pass
    # For each ligand / protein combination, we define the list 'combinations'
    # This list specifies all of the possible ligand/protein position 
    # combinations that are non-overlapping in a way that allows their
    # binding energy to be rapidly computed as follows:
    # 'combinations' is a list of the tuples '(pairlist, xshift, yshift, conf)'
    #  'xshift', 'yshift', and 'conf' specify the combination with
    #  the meaning described in the documentation string for this method.
    #  'pairlist' is a list of 2-tuples giving all pairs of ligand-protein
    #  residues that are in contact.  The 2-tuple '(i, j)' in 'pairlist'
    #  means that residue 'prot[i]' from the protein is in contact with 
    #  with residue 'ligand[j]' from the ligand.
    #
    # First, we look at the dictionary '_saved_combinations', which stores
    # the already saved combinations.
    # If we have more saved combinations than 'numsaved', get rid of half of them
    if numsaved < len(_saved_combinations):
        decorateditems = [(i[1][0], i) for i in _saved_combinations.iteritems()]
        decorateditems.sort()
        decorateditems.reverse()
        for i in decorateditems[int(numsaved / 2) : ]:
            del _saved_combinations[i[1][0]]
        for i in decorateditems[0 : int(numsaved / 2)]:
            _saved_combinations[i[1][0]] = (0, i[1][1][1])
    # If the combinations for these protein and ligand conformations 
    # are already stored, get them 
    savekey = (''.join(protconf), ''.join(ligandconf))
    try:
        _saved_combinations[savekey] = (_saved_combinations[savekey][0] + 1, _saved_combinations[savekey][1])
        combinations = _saved_combinations[savekey][1]
    except KeyError:
    # Combination is not stored, we will generate it and then store it.
        # First we create representations of the 2D lattices giving the
        # positions of all of the ligand residues, and all of the protein
        # residues.  'liganddict' is a dictionary which has a key at
        # '(x, y)' iff there is a ligand residue at position '(x, y)',
        # where the first ligand residue is placed at '(0, 0)'.  The
        # values for 'liganddict' are the residue numbers (0 <= number < length).
        # 'protdict' is the same for the protein residues.
        # We also compute 'minligandx', 'minligandy', 'maxligandx', 
        # and 'maxligandy' represeting the maximum and minimum x and
        # y coordinates for ligand residues.  Also the same for the
        # protein as 'minprotx', etc.
        dx = {'U':0, 'R':1, 'D':0, 'L':-1} # steps in x-direction
        dy = {'U':1, 'R':0, 'D':-1, 'L':0} # steps in y-direction
        rotate90conf = {'U':'L', 'R':'U', 'D':'R', 'L':'D'} # rotates 90 degrees
        #-------------------------------------------------------------
        def RotateLigand(in_conf):
            """Rotates ligand 90 degrees.

            Call is: '(conf, liganddict, minligandx, maxligandx, minligandy, 
                maxligandy) = RotateLigand(in_conf)'."""
            conf = ''.join([rotate90conf[c] for c in in_conf])
            ires = x = y = minligandx = maxligandx = minligandy = maxligandy = 0
            liganddict = {(x, y):ires}
            for c in conf:
                try:
                    x += dx[c]
                    y += dy[c]
                    ires += 1
                except KeyError:
                    raise ConformationsError, "Invalid ligand conformation in %r." % ligandconf
                liganddict[(x, y)] = ires
                minligandx = min(x, minligandx)
                maxligandx = max(x, maxligandx)
                minligandy = min(y, minligandy)
                maxligandy = max(y, maxligandy)
            if not (len(liganddict) == len(ligand) == ires + 1):
                raise ConformationsError, "Overlapping residues in ligand conformation %r." % ligandconf
            return (conf, liganddict, minligandx, maxligandx, minligandy, maxligandy)
        #-------------------------------------------------------------
        ires = x = y = minprotx = maxprotx = minproty = maxproty = 0
        protdict = {(x, y):0}
        for c in protconf:
            try:
                x += dx[c]
                y += dy[c]
                ires += 1
            except KeyError:
                raise ConformationsError, "Invalid protein conformation in %r." % protconf
            protdict[(x, y)] = ires
            minprotx = min(x, minprotx)
            maxprotx = max(x, maxprotx)
            minproty = min(y, minproty)
            maxproty = max(y, maxproty)
        if not (len(protdict) == len(prot) == ires + 1):
            raise ConformationsError, "Overlapping residues in protein conformation %r." % protconf
        # Now build all possible combinations.  Loop over all ligand rotations:
        combinations = []
        conf = ligandconf
        for rotation in range(4):
            (conf, liganddict, minligandx, maxligandx, minligandy, maxligandy) = RotateLigand(conf)            
            for xshift in range(minprotx - maxligandx - 1, maxprotx + 1 - minligandx):
                for yshift in range(minproty - maxligandy - 1, maxproty + 1 - minligandy):
                    pairlist = []
                    for ((x, y), jres) in liganddict.iteritems():
                        x += xshift
                        y += yshift
                        if (x, y) in protdict: # we have overlap
                            break # go to next position
                        for neighbor in ['U', 'D', 'R', 'L']:
                            try:
                                ires = protdict[(x + dx[neighbor], y + dy[neighbor])]
                                pairlist.append((ires, jres))
                            except KeyError:
                                pass
                    else: # no overlap, save this combination
                        combinations.append((pairlist, xshift, yshift, conf))
        _saved_combinations[savekey] = (0, combinations)
    # Loop over all of the combinations
    returntup = None
    for (pairlist, xshift, yshift, conf) in combinations:
        be = 0.0   # initialize binding energy to zero
        for (ires, jres) in pairlist:  # loop over interacting residues pairs
            try:
                be += interactions.miyazawa_jernigan["%s%s" % (prot[ires], ligand[jres])]
            except KeyError:
                raise ConformationsError, "Unrecognized residue in protein %r or ligand %r." % (prot, ligand)
        if returntup == None or be < returntup[0]:
            returntup = (be, xshift, yshift, conf)
    _saved_exactcombinations[exactsavekey] = (0, returntup)
    return returntup
#---------------------------------------------------------------------
def PrintConformation(seq, conf, file = sys.stdout, latex_format = False, ligand_tup = None, latex_justprot = False):
    """Prints the conformation to screen or file.

    Calls is: 'PrintConformation(seq, conf, [file = sys.stdout, 
        latex_format = False, ligand_tup = None])'
    'seq' is the protein sequence, as a string or list of 
        one letter amino acid codes.
    'conf' is the string specifying the conformation.  The string
        consists of the characters 'U', 'D', 'L', or 'R', and is of
        length 'len(seq)' - 1.  'conf[j]' is the conformation of bond  
        'j' (the bond between residues j and j + 1).  'U' means Up, 
        'D' means Down, 'R' means Right, and 'L' means Left.
    The conformation is printed to the output specified by 'file'.  By
        default, this is standard output: 'sys.stdout'.  If  
        'file' has a value different then it should be a writable file 
        object  In this case, the conformation
        is written to that file.  'file' must already be open for writing,
        and it is NOT closed by this method.
    'latex_format' is a boolean switch specifying if we should print the 
        conformation as a LaTex table.  If it is 'True', then we print
        it as a LaTex table, otherwise we do not.
    'latex_justprot' is a boolean switch specifying that if we
        are printing a LaTex table, we include no border or
        dots at empty lattice sites -- just the protein.
    'ligand_tup' can be used to specify that we also print a ligand to
        which the lattice protein is bound.  By default, 'ligand_tup'
        is 'None', meaning that there is no ligand bound to the 
        protein.  'ligand_tup' can also be set to the 4-tuple
        '(ligandseq, ligandconf, xshift, yshift)'.  'ligandseq'
        is a sequence giving the sequence of the ligand. 'ligandconf'
        is a sequence giving the conformation of the ligand in the
        same format as 'conf'.  'xshift' and 'yshift' specify where
        the ligand is located relative to the protein: move 'xshift'
        positions to the right and 'yshift' positions up from the
        position of the first protein residue to place the first
        ligand residue.  The ligand is printed in lower case letters."""
    # steps in directions for different conformations
    dx = {'U':0, 'R':2, 'D':0, 'L':-2}
    dy = {'U':2, 'R':0, 'D':-2, 'L':0}
    # Error check
    if len(seq) != len(conf) + 1:
        raise ConformationsError, "Lengths of 'seq' and 'conf' are incompatible."
    if len(seq) < 1:
        raise ConformationsError, "Sequence 'seq' is empty."
    # create a dictionary listing coordinates with residues at those sites
    # This dictionary is 'residue_coords'
    x = y = minx = maxx = miny = maxy = 0
    residue_coords = {(x, y):seq[0]}
    for ibond in range(len(conf)):
        newx = x + dx[conf[ibond]]
        newy = y + dy[conf[ibond]]
        # place the new coordinates
        residue_coords[(newx, newy)] = seq[ibond + 1]
        # place the bond
        if newx != x:
            residue_coords[((newx + x) / 2, y)] = '-'
        else:
            residue_coords[(x, (newy + y) / 2)] = '|'
        x = newx
        y = newy
        # update min/max coordinates
        minx = min(x, minx)
        miny = min(y, miny)
        maxx = max(x, maxx)
        maxy = max(y, maxy)
    # now place the ligand residues/bonds if we have a ligand
    if ligand_tup != None:
        # unpack the ligand tuple and do some error checking
        try:
            (ligandseq, ligandconf, xshift, yshift) = ligand_tup
        except (ValueError, TypeError):
            raise ConformationsError, "Invalid 'ligand_tup' of %r." % ligand_tup
        if len(ligandseq) != len(ligandconf) + 1:
            raise ConformationsError, "Ligand sequence of %r and ligand conformation of %r are of incompatible lengths." % (ligandseq, ligandconf)
        if not (isinstance(xshift, int) and isinstance(yshift, int)):
            raise ConformationsError, "Invalid 'xshift', 'yshift' pair of %r, %r." % (xshift, yshift)
        # make 'ligandseq' all lower case, as we will print lower case letters
        ligandseq = ''.join(ligandseq).lower()
        # now add the ligand residues/bonds to 'residue_coords'
        x = 2 * xshift
        y = 2 * yshift
        residue_coords[(x, y)] = ligandseq[0]
        # update the min/max coordinates
        minx = min(x, minx)
        miny = min(y, miny)
        maxx = max(x, maxx)
        maxy = max(y, maxy)
        for ibond in range(len(ligandconf)):
            # place the next residue
            newx = x + dx[ligandconf[ibond]]
            newy = y + dy[ligandconf[ibond]]
            residue_coords[(newx, newy)] = ligandseq[ibond + 1]
            # place the bond
            if newx != x:
                residue_coords[((newx + x) / 2, y)] = '-'
            else:
                residue_coords[(x, (newy + y) / 2)] = '|'
            x = newx
            y = newy
            # update the min/max coordinates
            minx = min(x, minx)
            miny = min(y, miny)
            maxx = max(x, maxx)
            maxy = max(y, maxy)
    # create the image
    xlength = maxx - minx + 5
    image = []
    if not latex_format: 
        # top buffer row
        for x in range(xlength):
            if x % 2:
                image.append(" ")
            else:
                image.append("*")
        image.append("\n")
        for x in range(xlength):
            image.append(" ")
        image.append("\n")
        # coordinates with structure
        for y in range(maxy, miny - 1, -1):
            if y % 2:
                image.append("  ")
                for x in range(minx, maxx + 1):
                    try:
                        image.append(residue_coords[(x, y)])
                    except KeyError:
                        image.append(" ")
                image.append("  ")
                image.append("\n")
            else:
                image.append("* ")
                for x in range(minx, maxx + 1):
                    if x % 2:
                        try:
                            image.append(residue_coords[(x, y)])
                        except KeyError:
                            image.append(" ")
                    else:
                        try:
                            image.append(residue_coords[(x, y)])
                        except KeyError:
                            image.append("*")
                image.append(" *")
                image.append("\n")
        # bottom buffer row
        for x in range(xlength):
            image.append(" ")
        image.append("\n")
        for x in range(xlength):
            if x % 2:
                image.append(" ")
            else:
                image.append("*")
        image.append("\n")
    else: # print in latex format
        if latex_justprot:
            emptysite = " &"
        else:
            emptysite = " $\\cdot$ &" # character for empty lattice sites
        vertbond = "$\\mid$ &" # vertical bonds
        horizbond = "--- &" # horizontal bonds
        emptybond = " &"  # empty bond
        endline = " \\\\\n" # end of line
        for (key, symbol) in residue_coords.iteritems():
            if symbol == '|':
                residue_coords[key] = vertbond
            elif symbol == '-':
                residue_coords[key] = horizbond
            elif len(symbol) == 1:
                residue_coords[key] = "%s &" % symbol
        # commands to create the table
        if latex_justprot:
            image.append("{\\begin{tabular}{")
        else:
            image.append("\\fbox{\\begin{tabular}{c")
        for i in range(xlength - 3):
            image.append("c@{}")
        image.append("ccc@{}}\n")
        # top buffer row
        if not latex_justprot:
            for x in range(xlength):
                if x % 2:
                    image.append(emptybond)
                else:
                    image.append(emptysite)
            image.append(endline)
            for x in range(xlength):
                image.append(emptybond)
            image.append(endline)
        # coordinates with structure
        for y in range(maxy, miny - 1, -1):
            if y % 2:
                if not latex_justprot:
                    image.append("%s%s" % (emptybond, emptybond))
                for x in range(minx, maxx + 1):
                    try:
                        image.append(residue_coords[(x, y)])
                    except KeyError:
                        image.append(emptybond)
                if not latex_justprot:
                    image.append("%s%s" % (emptybond, emptybond))
            else:
                if not latex_justprot:
                    image.append("%s%s" % (emptysite, emptybond))
                for x in range(minx, maxx + 1):
                    if x % 2:
                        try:
                            image.append(residue_coords[(x, y)])
                        except KeyError:
                            image.append(emptybond)
                    else:
                        try:
                            image.append(residue_coords[(x, y)])
                        except KeyError:
                            image.append(emptysite)
                if not latex_justprot:
                    image.append("%s%s" % (emptybond, emptysite))
            image.append(endline)
        # bottom buffer row
        if not latex_justprot:
            for x in range(xlength):
                image.append(emptybond)
            image.append(endline)
            for x in range(xlength):
                if x % 2:
                    image.append(emptybond)
                else:
                    image.append(emptysite)
            image.append(endline)
        # commands to end the table
        image.append("\\end{tabular}}\n")
    # make the image a string
    image = ''.join(image)
    file.write(image)
#---------------------------------------------------------------------
# End conformations.py
