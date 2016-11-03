from pynitefields import *
from balthasar.curve import Curve 
from balthasar.striations import Striations
import numpy as np
from functools import reduce
from operator import mul

class MUBs():
    """ Class to hold a complete set of mutually unbiased bases (MUBs)
        in a specified dimension. Currently, class MUBs fully supports only 
        qubit systems (i.e. systems of dimension :math:`2^n`).
        
        MUBs are commonly represented in two forms: collections of mutually
        unbiased basis vectors, or a table of disjoint, commuting operators 
        whose sets of mutual eigenvectors are mutually unbiased bases.
        We choose the latter representation, and represent MUBs as a table
        of tuples of the form (name, phase, matrix); the organization of this
        table will be specified after what follows.
        
        An operator in the MUB table is represented as a displacement 
        operator of the form 

        .. math::

         D(\\alpha, \\beta) = \Phi(\\alpha, \\beta) Z_\\alpha X_\\beta

        where :math:`\\alpha, \\beta` are field elements, and 
        :math:`\Phi(\\alpha, \\beta)` is a phase factor.
        The two operators Z and X are defined as having the action
        
        .. math::

              Z_\\alpha |\ell> = \omega(\\a \ell) |\ell>, 
              X_\\beta |\ell> = |\ell + \\beta>.
              
        For a single qubit these operators are exactly equal to the Pauli
        operators; for multiple qubits, tensor products thereof.
        These operators satisfy the Heisenberg-Weyl commutation relations, i.e.

        .. math::
              Z_\\alpha X_\\beta = \chi(\\alpha \\beta) X_\\beta Z_\\alpha

        where 

        .. math::
              
              \chi(\\alpha) = \exp(2 \pi i / p^n) ^ tr(\\alpha)

        is the character of the field, and the trace of a field element is 

        .. math::

            tr(\\alpha) = \\alpha + \\alpha^p + \cdots + \\alpha^{p^n - 1}

        in prime dimensions (where w represents the p^th root of unity), 
        and similarly in composite dimensions but w is replaced by the group
        character. For qubits Z and X are the simple Pauli operators. For
        qudits they are the generalized Paulis.
        These guys will be stored as tuples containing the phase factors
        separately because they're not always needed, i.e.
                      ( op_name,  phi(a, b), Z_a * X_b )
       
        By default, we will produce MUBs of the Desarguesian form, 
        where beta = lambda alpha.

        Attributes:
            field (GaloisField): The finite field of choice 
            p (int): A prime number, the dimension of a single particle
            n (int): The number of particles
            dim (int): The dimension of the system dim = :math:`p^n`
            w (pthRootOfUnity): The :math:`p^{th}` root of unity in the field.
            curves (list of Curves): The set of curves used to construct the table
            table: The table of operators, in form (name, phase, matrix)
            D (dictionary): Table of displacement operators, mapping of points
                            in phase space to the associated (phase, matrix)
            matrices (bool): Default true. Constructs the matrices associated
                             with the operators. If set to false, only the 
                             string representations of the operators is computed.

    """

    def __init__(self, f, **kwargs):
        """ Initialize a new table of MUBs.
            f is a GaloisField over which we should build the field.
              
            f must be represented in the self-dual basis, because it is
            this trace-orthogonality which allows us to properly factor the
            monomials in terms of single particle paulis.

            ----- Optional parameters in kwargs -----
            curves - The user can hand in a set of rays to generate MUBs from.
                     These must be a valid set of rays, which intersect only
                     at the point of origin (0, 0). If no curves are specified,
                     MUBs will be constructed using the Desarguesian rays.

            matrix - Tells whether or not to generate the matrix forms of the
                     operators or not. This is required to do things like 
                     generate/plot Wigner functions. This parameter is true
                     by default. However, for big systems, it takes a long 
                     time to generate them, and a user who is interested, say,
                     in only the operator table, or surviving coarse-grained 
                     operators, might not want to waste time generating them.
                     Set to 'false' to disable matrices.
        """
        if f.is_sdb() is False:
            print("Warning - you have passed a finite field whose \
                    expansion coefficients are not represented in \
                    the self-dual basis.")
            return None

        # ------------------------------------------------------------
        # Set some obvious parameters
        self.field = f # Keep a copy of the finite field to do math!
        self.p = f.p 
        self.n = f.n
        self.dim = f.dim
        self.w = f.w

        # Declaring these parameters here so we see them all together 
        self.curves = [] # Which curves to use
        self.table = [] # Operator table in operator form
        self.D = {} # Dictionary of displacement operators

        # Find and store the inverse of two only once, since we need
        # it in all our phases for odd primes
        self.twoinv = None
        if self.p != 2:
            if self.n == 1:
                self.twoinv = self.field[2].inv()
            else:
                for el in self.field:
                  if el.exp_coefs == ([2] + ([0] * (self.n -1))):
                    self.twoinv = el.inv()
        # ------------------------------------------------------------

        self.curves = []
        self.matrices = True

        # Set the curves, default to Desarguesian bundle if nothing passed in 
        if "curves" in kwargs: 
            if self.verify_curves(curves) == True:
                self.curves = curves 
        else: # Desarguesian curves
            self.curves = Striations.generate_rays(f)

        # The 
        if "matrix" in kwargs:
            if kwargs["matrix"] == False:
                self.matrices = False   

        # Build the operator table in matrix form at the same time
        self.table, self.D = self.build_operator_table()



    def phi(self, a, b):
        """ Phase function for the displacement operators.
            Handles every case (p odd/even, single/multi-qudit).
            Returns phase as either number, for qubits, or as a power
            of a pth root of unity for qudits.

            The phase condition in general for the displacement operators
            is that they must satisfy
                        phi^2(a, b) = gchar(-ab)
            Of course for qubits this becomes
                        phi^2(a, b) = gchar(ab)
            and ends up reducing to some powers of i.
        """
        if self.p == 2: # Qubits
            if self.n == 1: # Single qubit case, +/- i^ab
                return 1j ** (a * b).prim_power
            else: 
                return 1j ** tr(a * b)
        else: # Qudits
            if self.n == 1: # Single qudit case, w ^ (2^-1 ab)
                prefactor = (self.twoinv * a * b).prim_power
                return pow(self.w, prefactor) 
            else: # Multiple qudits
                return gchar(self.twoinv * a * b)


    def build_operator_table(self):
        """ Actually construct the operator table for the MUBs.
            This will be done using the displacement operators.

            We will consider our matrices of the form
                D(a, b) = phi(a, b) Z_a X_b 
            where a and b are field elements.

            We can further break this down by considering everything
            as an n-particle system, i.e.
                  b = sum (b_i theta_i), a = sum (a_i theta_i)
            for the expansions in the (almost) self-dual basis. Then,
                Z_a X_b = Z^a1 X^b1 \otimes ... \otimes Z^an X^bn
            where Z, X are the generalized Paulis.

            Note that this decomposition works only when there is 
            a true self-dual basis. For almost self-dual bases we will need
            to do something different.
        """
        table = []
        D = {}

        # Hold the generalized Paulis and the identity  
        Z = np.zeros((self.p, self.p), dtype=np.complex_)
        X = np.zeros((self.p, self.p), dtype=np.complex_)
        I = np.identity(self.p)
        
        if self.p == 2: # Simple case of qubit Paulis
            X = np.array([[0, 1], [1, 0]]) # X
            Z = np.array([[1, 0], [0, -1]]) # Z
        else:
            # Diagonal X, thanks to 
            # SO questions/10936767/rearranging-matrix-elements-with-numpy
            perm_order = [self.p - 1] + [x for x in range(self.p - 1)] 
            X = I[perm_order, :]

            # Diagonal generalized Z
            powers_of_w = [pow(self.w, i).eval() for i in range(self.dim)]
            np.fill_diagonal(Z, powers_of_w)

        # Now it's time to actually build the tuples of the operator table
        for curve in self.curves:
            row = [] # Each curve produces a different row of the table
            for point in curve: # (a, b)
                op= []
                op_mats = []

                a, b = point[0], point[1]

                if a == self.field[0] and b == self.field[0]:
                    continue # We ignore the identity

                phase = self.phi(a, b)
                z = a.exp_coefs # Expansion of the Z part
                x = b.exp_coefs # Expansion of the X part

                for idx in range(len(x)):
                    if z[idx] == 0 and x[idx] == 0: # Both coefs 0
                        op.append("I") # Tensor factor is identity
                    elif z[idx] == 0 and x[idx] != 0:
                        op.append("X" + ("" if x[idx] == 1 else str(x[idx])))
                    elif z[idx] != 0 and x[idx] == 0:
                        op.append("Z" + ("" if z[idx] == 1 else str(z[idx])))
                    else:
                        op.append("Z" + ("" if z[idx] == 1 else str(z[idx])) + \
                            "X" + ("" if x[idx] == 1 else str(x[idx])))

                    if self.matrices:
                        # Matrix for this chunk of the tensor product
                        Z_part = np.linalg.matrix_power(Z, z[idx])
                        X_part = np.linalg.matrix_power(X, x[idx])
                        op_mats.append(np.dot(Z_part, X_part))
                          
                if self.matrices:
                    # Tensor together all the matrices if the user wants it
                    matrix_op = reduce(np.kron, op_mats)
                    row.append( (op, phase, matrix_op) )
                
                    # Add the displacement operator to the matrix
                    D[(a, b)] = (phase, matrix_op)
                else:
                    row.append( (op, phase, None) )

            table.append(row) # Add rows to the table

            # Finally, add the identity operator to the D table if user wants
            if self.matrices:
                id_phase = self.phi(self.field[0], self.field[0])
                D[(self.field[0], self.field[0])] = (id_phase, np.identity(self.dim)) 

        return table, D


    def verify_curves(self, curves):
        """ TODO check that the properties of the provided curves are valid
        """
        if len(curves) != self.dim + 1:
            print("Error, not enough curves provided.")
            return False
        return True


    def print(self, matrix_form = False):
        np.set_printoptions(precision=4, threshold=np.nan, suppress=True)
        for row in self.table:
            for operator in row:
                print(str(operator[1]) + " ", end = "") # Phase
                print(" ".join(operator[0]) + "\t\t", end = "")
                if matrix_form: # Print as matrices
                    print()
                    print(operator[2]) # Matrix
            print("\n")
