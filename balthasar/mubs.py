from pynitefields import *
from balthasar.curve import Curve 
from balthasar.striations import Striations
import numpy as np
from functools import reduce

class MUBs():
    """ Class to hold a complete set of mutually unbiased bases
        in a specified dimension. 
        
        Member variables:
        p - The dimension of a single particle
        n - The number of particles
        dim - The dimension of the system p^n
        table - The table of operators
        matrix_table - The table of operators in explicit matrix form
        curves - The set of curves used to construct the table

        An operator in the MUB table is represented as a monomial of the
        form Z_alpha X_beta. By default, we will produce MUBs of the 
        Desarguesian form, where beta = lambda alpha.
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

        # Set some obvious parameters
        self.field = f # Keep a copy of the finite field to do math!
        self.p = f.p
        self.n = f.n
        self.dim = f.dim

        # If no curves have been passed in, we should just use the
        # Desarguesian bundle.

        self.curves = []
        if curves == []: # Default to Desarguesian curves
            striations = Striations(self.field)
            self.curves = striations.get_rays()
        else: # Curves specified by user
            if self.verify_curves(curves) == True:
                self.curves = curves 

        # Build the operator tables
        self.table = []
        self.matrix_table = []

        # Build the operator table in matrix form at the same time
        X = np.array([[0, 1], [1, 0]])
        Z = np.array([[1, 0], [0, -1]])
        ZX = - np.dot(Z, X)
        I = np.identity(2)
        matrix_dict = {"X" : X, "Z" : Z, "ZX" : ZX, "1" : I} # Easy access
    
        for curve in self.curves:
            row = []
            matrix_row = []
            for point in curve:
                op= []
                if point[0] == self.field[0] and point[1] == self.field[0]:
                    continue
                z = point[0].exp_coefs # Expansion of the Z part
                x = point[1].exp_coefs # Expansion of the X part
                for idx in range(len(z)):
                    if z[idx] == 0 and x[idx] == 0:
                        op.append("1")
                    elif z[idx] == 0 and x[idx] != 0:
                        op.append("X" + ("" if x[idx] == 1 else str(x[idx])))
                    elif z[idx] != 0 and x[idx] == 0:
                        op.append("Z" + ("" if z[idx] == 1 else str(z[idx])))
                    else:
                        op.append("Z" + ("" if z[idx] == 1 else str(z[idx])) + "X" + ("" if x[idx] == 1 else str(x[idx])))

                row.append(op)
                matrix_op = reduce(np.kron, (matrix_dict[i] for i in op)) # Compute matrix product
                matrix_row.append(matrix_op) 

            self.table.append(row) # Add to the tables
            self.matrix_table.append(matrix_row)

            

    def compute_generators(self):
        """ Compute the generators for each ray in the MUB table.
            Store two versions, one for the letter version and another for the matrices.
        """ 
        generators_operator = []
        generators_matrix = []
        
        for row_idx in range(len(self.table)):  # Compute for all rows of the table
            next_gen_op = self.table[row_idx][:2] # Collect first two elements
            next_gen_mat = self.matrix_table[row_idx][:2]
            
            num_gen = 2 # Current number of generators 
            next_idx = 2 # Index of next operator to check
            # Keep track of all possible products and add to this as we compute more generators
            generator_products = [next_gen_mat[0], next_gen_mat[1], np.dot(next_gen_mat[0], next_gen_mat[1])]

            while num_gen < self.field.n: # Keep going until we have n generators
                next_op = self.matrix_table[next_idx]
                
                # Check if the new matrix is in the span of the previous ones
                equality_test = [np.equal(next_op, gen_prod).all() for gen_prod in generator_products]

                if any(equality_test): # Invalid, increment and move on
                    next_idx += 1
                    continue 
                else: # If it's not found, it can be used as a generator, so update the span list
                    generator_products += [np.dot(next_op, gen_prod) for gen_prod in generator_products]
                    generator_products += [next_op]    

                    next_gen_op.append(self.table[row_idx][next_idx]) 
                    next_gen_mat.append(self.matrix_table[row_idx][next_idx]) 
                    num_gen += 1
                    next_idx += 1

            generators_operator.append(next_gen_op) 
            generators_matrix.append(next_gen_mat) 

        return generators_operator, generators_matrix


    def verify_curves(self, curves):
        # TODO check that the properties of the provided curves are valid
        if len(curves) != self.dim + 1:
            print("Error, not enough curves provided.")
            return False
        return True


    def print(self, matrix_form = False):
        if matrix_form: # Print as matrices
            np.set_printoptions(threshold=np.nan, suppress=True)
            for i in range(len(self.matrix_table)):
                for operator in self.table[i]:
                    print(" ".join(operator) + "\t\t", end = "")
                print("\n")
                for operator in self.matrix_table[i]:
                    print(operator) 
                    print("\n")
                print("\n")
        else: # Print as a table
            for row in self.table:
                for operator in row:
                    print(" ".join(operator) + "\t\t", end = "")
                print("\n")

