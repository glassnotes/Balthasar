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
        subfield_map - Map between elements of the large field with those
                       in the small field
        coarse_dim - A tuple representing the dimensions of the coarse system
        coset_reps - The coset representatives
        cosets - The cosets of the finite field (incl. representatives)
        coarse_D - The new (aggregate) displacement operators
        coarse_kernel - The operators at each point in the coarse phase space
    """

    def __init__(self, wf, coarse_field, **kwargs):
        """ Initialize a coarse-grained Wigner function based on a 
            fine-grained one.

            ----- Mandatory parameters -----
            wf - A fine-grained, i.e. normal, Wigner function
            coarse_field - The effective field after we've coarse-grained.
                           This allows us to specify how we want the 
                           coarse-graining to proceed.
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
        WignerFunction.__init__(self, wf.mubs)

        self.coarse_field = coarse_field

        # Loop through the optional arguments 
        self.basis = []
        self.cosets = []
        self.mode = "general"

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
                self.mode = kwargs["mode"]
                if self.mode != "general" and self.mode != "subfield":
                    print("Error, the mode must be 'general' or 'subfield'.")
                    return None
                # Dimension needs to be square for subfield mode. In other
                # words, the number of elements in each coset must be the same
                # as the size of each coset.
                if self.mode == "subfield":
                    if self.coarse_field.dim ** 2 != self.field.dim:
                        print("Error, dim must be square for subfield mode.")
                        return None

            if "cosets" in kwargs:
                # Make sure the cosets are valid
                field_as_list = [self.field[i] for i in range(self.dim)]
                cosets_flattened = [el for c in kwargs["cosets"] for el in c]
                if sorted(field_as_list) != sorted(cosets_flattened):
                    print("Error, are not a proper partitioning of the field.")
                    return None
                # If they're valid, choose representatives and carry on
                self.cosets = kwargs["cosets"]

        # Compute the subfield
        self.subfield, self.subfield_map = self.compute_subfield()
        if self.subfield == None:
            return None

        # If the cosets were not provided, compute them now. 
        if len(self.cosets) == 0:
            self.cosets = self.compute_cosets()

        # Compute the coarse displacement operators.
        self.coarse_D = self.compute_coarse_D()

        if self.mubs.matrices == False:
            print("Warning: coarse Wigner function constructed without matrices.")
            return

        # Using the coarse displacement operators, compute the coarse kernel
        # The points here are indexed by their coset representatives.
        self.coarse_kernel = self.compute_coarse_kernel()
        

    def compute_subfield(self):
        """ Compute the subfield from the big field. e.g. if our original
            system is dimension 16, and we are coarse-graining down to 
            dimension 4, we want to find the copy of F4 in F16.
        """
        subfield = []
        subfield_map = {}

        # For now let's concern ourselves only with the square case. 
        # Let the original field have dimension d. Then the subfield has
        # is generated by powers of the element f[sqrt(d) + 1]
        if self.coarse_field.dim ** 2 == self.field.dim:
            prim_element = self.field[self.coarse_field.dim + 1]
            subfield = [prim_element ** x for x in range(self.coarse_field.dim)]    

            # Make the subfield map
            for i in range(len(subfield)):
                subfield_map[subfield[i]] =  self.coarse_field[i]
        else:
            # If the dimension is not square, we can only pull out a copy of
            # what Luis affectionately calls the "mother field", i.e. the
            # base prime field. These are the elements whos exp_coefs are all
            # 0 except for the first one. Find these!

            # We need the field in the polynomial basis to do this properly,
            # because in the self-dual basis the elements are not right form.
            field_poly_basis = GaloisField(self.field.p, self.field.n, \
                self.field.coefs)

            for el in field_poly_basis: 
                if all(a == 0 for a in el.exp_coefs[1:]):
                    el_in_sdb = self.field[el.prim_power] # "Real" element
                    subfield.append(el_in_sdb)
                    subfield_map[el_in_sdb] = self.coarse_field[el.exp_coefs[0]]

        return subfield, subfield_map


    def compute_polynomial_basis(self):
        """ Return the polynomial basis for coarse graining. 
            When the dimension isn't square, this 
            is just the polynomial basis in the 'big' field. When the dim     
            is square we choose the polynomial basis of the big field wrt     
            the small field. For example, in dimension 16 coarse-graining to  
            dimension 4, we choose the basis {1, \sigma} rather than the      
            full polynomial basis because dimension 16 is 2-dim vector space  
            over the dim 4 case."""
        if self.coarse_field.dim ** 2 == self.field.dim: # Square dimensions
            return [self.field[-1], self.field[1]]
        else: # Standard polynomial basis
            return [self.field[-1]] + [self.field[x] for x in range(1, self.field.n)] 


    def compute_cosets(self):
        """ Compute a set of coset representatives if the user doesn't 
            provide one. Two cases to consider here, one where there is
            a basis provided to make the cosets, the other where we use the
            subfield.
        """
        cosets = []

        if self.mode == "general":
            # In the general mode, for a system p^mn, we can really only take
            # out one copy of p in the most general case by using this method.
            # However, we can also apply the general mode to square dimensions
            # By choosing a 'small' basis and use the subfield elements
            # as coset representatives.
            first_coset = []

            # First generate all combinations of the coefficients with the
            # subfield except the first element (which is 1 in most cases).
            coefficient_lists = product(self.subfield, repeat = len(self.basis) - 1) 
            for c in coefficient_lists:
                l_comb = [c[i-1] * self.basis[i] for i in range(1, len(c) + 1)]
                s = self.field[0]
                for el in l_comb:
                    s += el
                first_coset.append(s)

            for i in range(len(self.subfield)):
                next_coset = [self.subfield[i] * self.basis[0]
                    + first_coset[t] for t in range(len(first_coset))]
                cosets.append(next_coset)

        else:
            l_om = self.field[1] # Primitive element of original field
            b_om = self.subfield[1] # Primitive element of subfield

            for c_idx in range(len(self.subfield)):
                if c_idx == 0:
                    cosets.append(self.subfield)
                elif c_idx == 1: # Special case because of how we defined pow
                    next_coset = [x + l_om for x in self.subfield]
                    cosets.append(next_coset)
                else:
                    next_coset = [x + l_om * (b_om ** (c_idx - 1)) \
                        for x in self.subfield]
                    cosets.append(next_coset)

        return cosets


    def compute_coarse_D(self):
        """ Compute the coarse-grained displacement operators. 
            This will be a subset of the fine-grained displacement operators
            which 'survive' a sum.
        """
        coarse_D = {}

        if self.field.p == 2:
            survivors = []

            for alpha in self.field:
                l = [(-1) ** tr(self.cosets[0][i] * alpha) \
                    for i in range(len(self.cosets[0]))]
                if sum(l) != 0:
                    survivors.append(alpha.prim_power)
            print(survivors)

            # Collect the surviving operators into a table
            # Note that the MUB table does not contain identity operators
            # at the beginning of each row - thus, we will have to take only 
            # the non-zero survivors, and change the index by -1.
            surviving_ops = {}
            for slope in self.subfield:
                table_row = self.mubs.table[slope.prim_power]
                surviving_ops[slope.prim_power] = \
                    [table_row[x-1][0] for x in survivors[1:]]

            # Infinite slope case
            infty_row = self.mubs.table[-1]
            surviving_ops["inf"] = [infty_row[x-1][0] for x in survivors[1:]]

        else:
            print("Sorry, qudit coarse-graining not currently implemented.")
            return None

        return surviving_ops 


    def compute_coarse_kernel(self):
        """ Compute the kernel of the coarse-grained Wigner function. 
            This is done by 'globbing' together the fine-grained kernel coset
            by cosets. The final kernel will be indexed by points in the 
            subfield, though really I suppose it should be indexed by 
            coset representatives in a way.

            Returns a dictionary mapping a pair of coset reps (a, b) to W(a,b).
        """

        ckernel = {}

        # Coarse kernel should be indexed by cosets / subfield elements
        for alpha in range(len(self.subfield)):
            for beta in range(len(self.subfield)):
                coarse_point = (self.coarse_field[alpha], \
                    self.coarse_field[beta])

                mat = np.zeros((self.field.dim, self.field.dim), \
                    dtype = np.complex_)
                
                # Sum up the kernel points from the appropriate cosets
                for x in self.cosets[alpha]:
                    for y in self.cosets[beta]:
                        fine_point = (x, y)
                        mat = mat + self.kernel[fine_point]

                ckernel[coarse_point] = (1.0 / self.coarse_field.dim) * mat
                
        return ckernel
                


    def compute_wf(self, state):                                                
        """ Compute the probabilities in the coarse Wigner function for a 
             given state. 
             Input: state, a numpy array representing either a ket or a 
             density matrix.
        """                                                                     
        # For the coarse Wigner function the dimension is that of the 
        # underlying affine plane.
        W = np.zeros((self.coarse_field.dim, self.coarse_field.dim)) 

        # Turn kets into density operators if need be.
        if state.shape[0] == 1:                                                 
            state = np.outer(state, np.conj(state))                             

        # The coarse Wigner function is indexed by the subfield, so use this.
        for a in range(self.coarse_field.dim):
            for b in range(self.coarse_field.dim):
                a_in_cf = self.subfield_map[self.subfield[a]] 
                b_in_cf = self.subfield_map[self.subfield[b]] 
                index = (a_in_cf, b_in_cf)
                mat = np.trace(np.dot(state, self.coarse_kernel[index]))              
                W[a][b] = (1.0 / self.coarse_field.dim) * mat          

        return W   
