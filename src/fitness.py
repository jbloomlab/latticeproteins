#!/usr/bin/python
# Begin fitness.py
#---------------------------------------------------------------------
"""Module for calculating fitnesses of lattice protein sequences.

Written by Jesse Bloom, 2004."""
#----------------------------------------------------------------------
import math, conformations, sys
#----------------------------------------------------------------------
class FitnessError(Exception):
    """Error computing lattice protein fitness."""
#----------------------------------------------------------------------
class Fitness(object):
    """Instances of this class compute fitnesses of lattice protein sequences."""
    #------------------------------------------------------------------
    def __init__(self, temp, conformations, dGdependence, targets, ligand = None, nofitness = -1.0e10):
        """Creates a new instance of fitness evaluator.

        Call is: 'f = Fitness(temp, conformations, dGdependence, targets,
            [ligand = None, nofitness = -1e10])'
        'temp' is the temperature at which the fitness is computed.
        'conformations' is the 'conformations.Conformations' object
            used to fold the protein sequences.  'conformations.Length()'
            specifies the length of the protein sequences that can be
            folded.
        'dGdependence' specifies how the fitness depends on the free
            energy of folding of the sequence:
            * if 'dGdependence' is number, then any sequence with a
                free energy of folding <= 'dGdependence' has a fitness
                of one, and any sequence with a free energy of folding
                > 'dGdependence' has a fitness of 'nofitness'. 
            * if 'dGdependence' is the string 'fracfolded' then the fitness
                of the sequence is the fraction of the sequences that
                will be folded at temperature 'temp' at equilibrium.
                If 'dG' is the free energy of folding, this fraction is:
                'f = 1 / (1 + exp(dG / temp)'
            * if 'dGdependence' is the string 'negstability', then the
                fitness of the sequence is just negative one times the
                stability of the sequence.  The negative one is so
                that more stable sequences have higher fitnesses.
        'targets' specifies information about the target conformation(s)
            to which the sequence is folded as follows:
            * If 'targets' is 'None' then the protein is folded to its
                lowest energy conformation, whatever this is.
            * If 'targets' is a string specifying a specific conformation
                then the free energies of folding are for folding the 
                protein to this specific conformation.
            * If 'targets' is an integer, then this integer specifies the
                number of contacts which the lowest energy conformation
                must have.  If the lowest energy conformation has a
                number of contacts different from 'targets', the
                fitness is 'nofitness'.
        'ligand' is an optional argument that is used if we are looking
            for a protein that binds a ligand.  By default, it is 'None'
            meaning that no ligand binding is considered.  If it is set 
            to another value, it should be the 3-tuple '(ligand, 
            ligandconf, stabcutoff)' where 'ligand' and 'ligandconf'
            are both strings describing a ligand as detailed in the
            documentation string for the 'conformations.BindLigand'
            method.  'stabcutoff' is a number specifying the
            stability cutoff for the protein to fold.  In this case,
            'dGdependence' no longer has any meaning.  The protein
            is folded according to the parameters 'temp' and 
            'targets' as normal.  If the free energy of folding of
            the protein is > 'stabcutoff', then the returned fitness is 
            zero.  If the free energy of folding of the protein is <=
            'stabcutoff', then the returned fitness 'exp(-be)' where
            'be' is the binding energy of the ligand to the protein
            in the folded conformation.
        'nofitness' is an optional argument specifying the fitness returned 
            for a sequence that does not satisify the requirements set
            by 'dGdependence' or 'targets'.  By default, it is the very
            negative number -1.0e10."""
        # Assign class instance variables and error check
        self._temp = temp
        if not (isinstance(self._temp, (int, float)) and temp > 0):
            raise FitnessError("Invalid 'temp' of %r." % temp)
        self._conformations = conformations
        self._dGdependence = dGdependence
        if not (isinstance(dGdependence, (int, float)) or dGdependence == 'fracfolded' or dGdependence == 'negstability'):
            raise FitnessError("Invalid 'dGdependence' of %r." % dGdependence)
        self._targets = targets
        if not (targets == None or (isinstance(targets, str) and len(targets) == self._conformations.Length() - 1) or isinstance(targets, int)):
            raise FitnessError("Invalid 'targets' of %r." % targets)
        self._ligand = ligand
        if ligand == None:
            pass
        elif isinstance(ligand, tuple) and len(ligand) == 3:
            (lig, ligconf, stabcutoff) = ligand
            if not (isinstance(lig, str) and isinstance(ligconf, str) and len(lig) == len(ligconf) + 1):
                raise FitnessError("%r does not specify a valid ligand." % ligand)
            if not (isinstance(stabcutoff, (int, float))):
                raise FitnessError( "Invalid 'stabcutoff' of %r." % stabcutoff)
        else:
            raise FitnessError("Invalid 'ligand' of %r." % ligand)
        self._nofitness = nofitness
    #---------------------------------------------------------------------
    def Fitness(self, seq):
        """Computes the fitness of a sequence.

        Call is: 'x = f.Fitness(seq)'
        'seq' is a list or string specifying a protein sequence of length
            'f.Length()'.
        'x' is returned as the fitness of this protein sequence."""
        if len(seq) != self.Length():
            raise FitnessError("Invalid 'seq' of %r." % seq)
        if isinstance(self._targets, str):
            # folding to a target conformation
            (dG, conf, numcontacts) = self._conformations.FoldSequence(seq, self._temp, self._targets) 
        else:
            # folding to lowest energy conformation
            (dG, conf, numcontacts) = self._conformations.FoldSequence(seq, self._temp) 
        if isinstance(self._targets, int) and self._targets != numcontacts:
            # wrong number of contacts
            return self._nofitness
        if self._ligand:
            if dG > self._ligand[2]:
                return 0 # does not stably fold
            else:
                be = conformations.BindLigand(seq, conf, self._ligand[0], self._ligand[1])[0]
                return math.exp(-be) # compute the fitness
        elif self._dGdependence == 'fracfolded':
            f = 1.0 / (1.0 + math.exp(dG / self._temp))
            return f
        elif self._dGdependence == 'negstability':
            return -dG
        else:
            # free energy cutoff
            if dG <= self._dGdependence:
                return 1.0
            else:
                return self._nofitness 
    #---------------------------------------------------------------------
    def Stability(self, seq):
        """Computes the stability of a sequence if it is below cutoff.

        Call is: 'dGf = f.Stability(seq)'
        'seq' is the sequence we are folding.  
        If 'dGdependence' is set to a free energy cutoff, then if 
            the dGf > dGdependence, dGf is just returned as 'None'.
        'dGf' is the free energy of folding of the sequence to the 
            target conformation."""
        if isinstance(self._dGdependence, (int, float)):
            return self._conformations.FoldSequence(seq, self._temp, self._targets, dGf_cutoff = self._dGdependence)[0]
        else:
            return self._conformations.FoldSequence(seq, self._temp, self._targets)[0]
    #---------------------------------------------------------------------
    def Info(self, file = sys.stdout):
        """Prints information about the fitness function.

        Call is: 'f.Info([file = sys.stdout])'
        A summary of the fitness function is printed.
        'file' is an optional argument specifying where the information
            is printed.  By default, it is standard output ('sys.stdout').
            If it is set to another value, it must be an open file-like
            object.  'file' is NOT closed by this method after writing."""
        file.write("Fitnesses are computed for proteins of length %d at a temperature of %.3f.\n" % (self.Length(), self._temp))
        if self._ligand:
            file.write("The fitness of a protein is equal to exp(-be) where be is the binding energy of the protein to ligand %s in conformation %s if the protein folds with a free energy of folding <= %4.f, and zero otherwise." % self._ligand)
        elif self._dGdependence == 'fracfolded':
            file.write("The fitness of a protein is equal to the fraction of the proteins that are folded.\n")
        elif self._dGdependence == 'negstability':
            file.write("Fitness is negative one times the free energy of folding.\n")
        else:
            file.write("A protein has a fitness of %r if its free energy is greater than %r, and a fitness of one otherwise.\n" % (self._nofitness, self._dGdependence))
        if isinstance(self._targets, str):
            file.write("Free energies of folding are to the target conformation %r.\n" % self._targets)
        elif isinstance(self._targets, int):
            file.write("Fitnesses are calculated as described above only for proteins with lowest energy conformation having %d contacts.  Otherwise the fitness is %r.\n" % (self._targets, self._nofitness))
        else:
            file.write("The free energy of folding is to the lowest energy conformation.\n")
    #---------------------------------------------------------------------
    def Length(self):
        """Returns the sequence length for which fitnesses are computed."""
        return self._conformations.Length()
    #---------------------------------------------------------------------
#---------------------------------------------------------------------------
# End fitness.py
