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

    def __init__(self, wf, dimensions, **kwargs):
        """ Initialize a coarse-grained Wigner function based on a 
            fine-grained one.

            ----- Mandatory parameters -----
            wf - A fine-grained, i.e. normal, Wigner function
            coarse_dims - The dimensions over which to partition the
                          field, e.g. for coarse graining dimension 16
                          you would pass (4, 4), for 8 (2, 4), etc.
            --------------------------------

            ----- Optional parameters in kwargs -----
            basis - Specify an alternate basis to do coarse graining over.
                    Default is the polynomial basis.
            cosets - A set of cosets to use. Will default to a standard
                     choice unless user provides a full set.
            mode - 'general' or 'subfield'. General mode works for any dim,
                    but subfield mode works only for square dimensions. 
            --------------------------------
        """
        # Initialize all the standard parameters
        WignerFunction.__init__(wf.mubs)

        # The dimensions of the system are as follows. The grid will be
        # partitioned into a len(cosets) x len(cosets) grid; each chunk
        # of that grid will have dimensions equal to the number of
        # elements within each coset.
        # The effective affine space is of dimension self.coarse_dims[0]
        self.coarse_dims = dimensions
        if self.coarse_dims[0] * self.coarse_dims[1] != self.dim:
            print("Error, invalid coarse grained space dimensions.")
            print("Dimensions must multiply to the dimension of the system.")
            return None

        # Compute the subfield F_dimensions[0] in the big field

        # Loop through the optional arguments 
        self.basis = []
        self.coset_reps = []
        self.cosets = []
        self.mode = 'general'

        # We will probably need to do some better validity checks here but this
        # should be okay for now (hopefully...)
        if kwargs is not None:
            # Choose the basis with which to do coarse-graining (usual this 
            # will matter only in the general case.
            if "basis" in kwargs:
                self.basis = kwargs["basis"]
            else:
                self.basis = self.compute_polynomial_basis()

            if "mode" in kwargs:
                self.mode = kwargs[mode]
                if self.mode != "general" and self.mode != "subfield":
                    print("Error, the mode must be 'general' or 'subfield.")
                    return None
                # Dimension needs to be square for subfield mode. In other
                # words, the number of elements in each coset must be the same
                # as the size of each coset.
                if self.mode == "subfield":
                    if self.dimensions[0] != self.dimensions[1]:
                        print("Error, dim must be square for subfield mode.")
                        return None

            if "cosets" in kwargs:
                # Make sure the cosets are valid
                field_as_list = [f[i] for i in range(self.dim)]
                cosets_flattened = [el for el in c for c in kwargs["cosets"]]
                if sorted(field_as_list) != sorted(cosets_flattened):
                    print("Error, are not a proper partitioning of the field.")
                    return None
                # If they're valid, choose representatives and carry on
                self.coset_reps = [c[0] for c in cosets]
                self.cosets = cosets

        # If the cosets were not provided, compute them now. 
        if len(self.cosets) == 0:
            self.coset_reps, self.cosets = self.compute_cosets()

        # Compute the coarse displacement operators.
        self.coarse_D = self.compute_coarse_displacement_operators()

        # Using the coarse displacement operators, compute the coarse kernel
        # The points here are indexed by their coset representatives.
        self.coarse_kernel = self.compute_coarse_kernel()
        

    def compute_polynomial_basis(self):
        # Return the polynomial basis
        return [self.field[x] for x in range(1, self.field.n)] + [self.field[-1]]


    def compute_cosets(self):
        """ Compute a set of coset representatives if the user doesn't 
            provide one. Two cases to consider here, one where there is
            a basis provided to make the cosets, the other where we use the
            subfield.
        """

        if self.mode == "general":
            # In the general mode, for a system p^mn, we can really only take
            # out one copy of p in the most general case by using this method.
            first_coset = []

            # Compute all linear combinations of all but 1 using the subfield 
            # First generate all combinations of the coefficients
            coefficient_lists = product(self.subfield, repeat = self.field.n - 1) 
            for coeffs in coefficient_list:
                l_comb = [coeffs[i] * self.basis[i] for i in range(len(coeffs))]
                s = self.field[0]
                for el in l_comb:
                    s += el
                first_coset.append(s)

            # Use the remaining element 

            return first_coset
        else:
            return None

    def compute_coarse_D(self):




    def compute_coarse_kernel(self):
        """ Compute the kernel of the coarse-grained Wigner function. This
            function is pretty much the same as the one for the normal WF, 
            but we need to index the kernel by coset representatives rather
            than plain old field elements this time.

            Returns a dictionary mapping a pair of coset reps (a, b) to w(a,b).
        """

        ckernel = {}

        ckernel_00 = np.zeros((self.dim, self.dim), dtype = np.complex_)
        if self.field.p != 2:
            for key in self.coarse_D.keys():
                ckernel_00 = ckernel_00 + \
                    (self.coarse_D[key][0].eval() * self.coarse_D[key][1]
            ckernel_00 = ckernel_00 / self.dim
            ckernel[(self.coset_reps[0], self.coset_reps[0])] = ckernel_00
        else:
            if self.field.n == 1:
                for key in self.coarse_D.keys():
                    ckernel_00 = ckernel_00 + \
                        (self.coarse_D[key][0] * self.coarse_D[key][1])
                ckernel_00 = ckernel_00 /self.dim
                ckernel[(self.coset_reps[0], self.coset_reps[0])] = ckernel_00

        return ckernel
                


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
