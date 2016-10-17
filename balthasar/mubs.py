from pynitefields import *
from balthasar.curve import Curve 
from balthasar.striations import Striations
import numpy as np
from functools import reduce
from operator import mul

class MUBs():
    """ Class to hold a complete set of mutually unbiased bases
        in a specified dimension. 
        
        Member variables:
        _field - The finite field of choice 
        _p - The dimension of a single particle
        _n - The number of particles
        _d - The dimension of the system _d = _p^_n
        _table - The table of operators
        _matrix_table - The table of operators in explicit matrix form
        _curves - The set of curves used to construct the table

        An operator in the MUB table is represented as a displacement 
        operator of the form 
                          D(a, b) = phi(a, b) U_b V_a
        where a, b are field elements, and phi(a, b) is a phase factor.
        The two operators U and V are defined as 
                  U^a |l> = |l + a>,   V^b |l> = w(bl) |l>
        in prime dimensions (where w represents the p^th root of unity), 
        and similarly in composite dimensions but w is replaced by the group
        character. I'm calling them U and V because in general they are the
        elements of the Heisenberg-Weyl group and not just the qubit Pauli X/Z.
        These guys will be stored as tuples containing the phase factors
        separately because they're not always needed, i.e.
                           ( op_name,  phi(a, b), U_b * V_a )
       
        By default, we will produce MUBs of the Desarguesian form, 
        where beta = lambda alpha.
    """

    """ The following series of functions provides phase factors for
        the displacement operators depending on what case we're in.
    """


    def __init__(self, f, curves = []):
        """ Initialize a new table of MUBs.
            f is a GaloisField over which we should build the field.
            Ideally, f is represented in self-dual basis form.
        """
        if f.is_sdb() is False:
            print("Warning - you have passed a finite field whose \
                    expansion coefficients are not represented in \
                    the self-dual basis.")

        # ------------------------------------------------------------
        # Set some obvious parameters
        self.field = f # Keep a copy of the finite field to do math!
        self.p = f.p 
        self.n = f.n
        self.d = f.dim
        self.w = f.w

        # Declaring these parameters here so we see them all together 
        self.curves = [] # Which curves to use
        self.table = [] # Operator table in operator form
        self.D = {} # Dictionary of displacement operators
        # ------------------------------------------------------------

        # Set the curves, default to Desarguesian bundle if nothing passed in 
        if curves == []: 
            striations = Striations(self.field)
            self.curves = striations.get_rays()
        else: # Curves specified by user
            if self.verify_curves(curves) == True:
                self.curves = curves 


        # Build the operator table in matrix form at the same time
        self.table, self.D = self.build_operator_table()



    def phi(self, a, b):
        """ Phase function for the displacement operators.
            Handles every case (p odd/even, single/multi-qudit).
            Returns phase as either number, or pth root of unity.
        """
        if self.p == 2: # Qubits
            if self.n == 1: # Single qubit case, +/- i^ab
                return 1j ** (a * b).prim_power
            else:
                return 1 # TODO
        else: # Qudits
            if self.n == 1: # Single qudit case, w ^ (2^-1 ab)
                prefactor = (self.field[2].inv() * a * b).prim_power
                return pow(self.w, prefactor) 
            else: # Multiple qudits
                two = None # First find two and it's inverse
                for el in self._field:
                    if el.exp_coefs == ([2] + ([0] * self._field.n - 1)):
                        two = el
                          
                if two == None:
                    print("Welp, something went wrong.")
                    return None

                return gchar(two.inv() * a * b)


    def build_operator_table(self):
        """ Actually construct the operator table for the MUBs.
            This will be done using the displacement operators.

            We will consider our matrices of the form
                D(a, b) = phi(a, b) U_b V_a
            where a and b are field elements.

            We can further break this down by considering everything
            as an n-particle system, i.e.
                  b = sum (b_i theta_i), a = sum (a_i theta_i)
            for the expansions in the (almost) self-dual basis. Then,
                U_b V_a = U^b1 V^a1 \otimes ... \otimes U^bn V^bn
            where U, V are the generalized Paulis (shift/diagonal).
            There are also some coefficients in here we'll have to figure out.
        """
        table = []
        D = {}

        # Hold the generalized Paulis and the identity  
        U = np.zeros((self.p, self.p), dtype=np.complex_)
        V = np.zeros((self.p, self.p), dtype=np.complex_)
        I = np.identity(self.p)
        
        if self.p == 2: # Simple case of qubit Paulis
            U = np.array([[0, 1], [1, 0]]) # X
            V = np.array([[1, 0], [0, -1]]) # Z
        else:
            # Diagonal X, thanks to 
            # SO questions/10936767/rearranging-matrix-elements-with-numpy
            perm_order = [self.p - 1] + [x for x in range(self.p - 1)] 
            U = I[perm_order, :]

            # Diagonal generalized Z
            powers_of_w = [pow(self.w, i).eval() for i in range(self.d)]
            np.fill_diagonal(V, powers_of_w)

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
                u = b.exp_coefs # Expansion of the U part
                v = a.exp_coefs # Expansion of the V part

                for idx in range(len(v)):
                    if u[idx] == 0 and v[idx] == 0: # Both coefs 0
                        op.append("I") # Tensor factor is identity
                    elif u[idx] == 0 and v[idx] != 0:
                        op.append("V" + ("" if v[idx] == 1 else str(v[idx])))
                    elif u[idx] != 0 and v[idx] == 0:
                        op.append("U" + ("" if u[idx] == 1 else str(u[idx])))
                    else:
                        op.append("V" + ("" if v[idx] == 1 else str(v[idx])) + \
                            "U" + ("" if u[idx] == 1 else str(u[idx])))

                    # Matrix for this chunk of the tensor product
                    V_part = np.linalg.matrix_power(V, v[idx])
                    U_part = np.linalg.matrix_power(U, u[idx])
                    op_mats.append(np.dot(V_part, U_part))
                        
                # Tensor together all the matrices 
                matrix_op = reduce(np.kron, op_mats)
                
                # Append the tuple to the row
                row.append( (op, phase, matrix_op) )

                # Add the displacement operator to the matrix
                D[(a, b)] = (phase, matrix_op)

            table.append(row) # Add to the tables

            # Finally, add the identity operator to the D table
            D[(self.field[0], self.field[0])] = (self.w ** 0, np.identity(self.d)) 

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
