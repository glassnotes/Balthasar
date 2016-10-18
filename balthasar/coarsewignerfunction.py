from pynitefields import *
from balthasar.wignerfunction import WignerFunction
import numpy as np

from itertools import product 

class CoarseWignerFunction(WignerFunction):
    """ Class to store and plot a coarse-grained discrete Wigner function.
        These are generated from fine-grained discrete Wigner functions by
        choosing some sort of coset partitioning of the underlying field.

        CoarseWignerFunction inherits from WignerFunction, since once the
        construction is complete, a lot of the functionality remains the same,
        e.g. plotting, computing the marginals.

        Parameters (inherited):
        ===========
        field - The finite field over which the Wigner function is defined
        dim - The dimension of the system
        mubs - The MUBs associated with this Wigner function
        D- The dictionary of displacement operators
        kernel - Operators at each point in discrete phase space
        
        Parameters (new):
        ===========
        subfield - The effective subfield within the large field
        coarse_dim - A tuple representing the dimensions of the coarse system
        coset_reps - The coset representatives
        cosets - The cosets of the finite field (incl. representatives)
        coarse_D - The new (aggregate) displacement operators
        coarse_kernel - The operators at each point in the coarse phase space
    """

    def __init__(self, wf, dimensions, cosets = []):
        """ Initialize a coarse-grained Wigner function based on a 
            fine-grained one.

            wf - A fine-grained, i.e. normal, Wigner function
            cosets - A set of cosets to use. Will default to a standard
                     choice unless user provide a full set.
                     TODO implement something that allows choice between
                          normal/special case.
        """
        # Initialize all the standard parameters
        WignerFunction.__init__(wf.mubs)

        # Compute the subfield F_dimensions[0] in the big field


        # Initialize the coset reps and compute the rest of the cosets
        self.cosets = cosets

        if len(self.cosets) == 0:
            self.coset_reps, self.cosets = self.compute_cosets()
        else:
            self.coset_reps = self.cosets[0]

        # The dimensions of the system are as follows. The grid will be
        # partitioned into a len(cosets) x len(cosets) grid; each chunk
        # of that grid will have dimensions equal to the number of
        # elements within each coset.
        # The effective affine space is of dimension self.coarse_dim[0]
        self.coarse_dim = (dimensions[0], dimensions[1])

        # Compute the coarse displacement operators.

        # Using the coarse displacement operators, compute the coarse kernel
        # The points here are indexed by their coset representatives.

    def compute_cosets(self):
        """ Compute a set of coset representatives if the user doesn't 
            provide one.
        """
        first_coset = []

        # Start with the polynomial basis
        p_basis = [self.field[x] for x in range(1, self.field.n)] + \
                     [self.field[-1]] # Make sure to add element 1
    
        # Compute all linear combinations of all but 1 using the subfield 
        # First generate all combinations of the coefficients
        coefficient_lists = product(self.subfield, repeat = self.field.n - 1) 
        for coeffs in coefficient_list:
            linear_comb = [coeffs[i] * p_basis[i] for i in range(len(coeffs))]
            s = 0
            for el in linear_comb:
                s += el
            first_coset.append(s)

        # Use the remaining elements

        return first_coset

    def compute_cosets(self):
        """ Compute the cosets based on the coset representatives.
        """
        # The number of coset reps should divide the dimension perfectly.
        # Of course there is more to it than this, but we can handle it later
        if self.dim % len(self.coset_reps) != 0:
            print("Error, invalid coset representatives.")
            return 

        



    def compute_wf(self, state):                                                
        """ Compute the probabilities in the coarse Wigner function for a 
             given state. 
             Input: state, a numpy array representing either a ket or a 
             density matrix.
        """                                                                     
        # For the coarse Wigner function the dimension is that of the 
        # underlying affine plane.
        W = np.zeros((self.coarse_dim[0], self.coarse_dim[0])) 

        # Turn kets into density operators if need be.
        if state.shape[0] == 1:                                                 
            state = np.outer(state, np.conj(state))                             

        # The coarse Wigner function is indexed by cosets / coset reps, so
        # loop over these to compute stuff.
        for a in range(len(self.coset_reps)):
            for b in range(len(self.coset_reps)):
                index = (self.coset_reps[a], self.coset_reps[b])
                mat = np.trace(np.dot(state, self.kernel[index]))              
                W[a][b] = (1.0 / self.dim) * mat          

        return W   
