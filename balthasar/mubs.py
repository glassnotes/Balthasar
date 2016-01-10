from pynitefields import *

class MUBs():
    """ Class to hold a complete set of mutually unbiased bases
        in a specified dimension. 
        
        Member variables:
        p - The dimension of a single particle
        n - The number of particles
        dim - The dimension of the system p^n
        table - The table of operators
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
        if curves == []:
            self.curves = [[0, self.field[i]] for i in range(self.dim)] 
        else:
            self.curves = curves 

        # Build the operator table
        self.table = []

    
        for curve_idx in range(len(self.curves) + 1):
            row = []
            if curve_idx == len(self.curves): # Append the horizontal rows (alpha = 0)
                ordered_pairs = [(self.field[0], self.field[i]) for i in range(1, self.field.dim)]
            else: # Evaluate the function - don't need to consider the (0, 0) case
                ordered_pairs = [(el, self.field.evaluate(self.curves[curve_idx], el)) for el in self.field][1:]

            for pair in ordered_pairs:
                operator = []
                z_coords = pair[0].exp_coefs
                x_coords = pair[1].exp_coefs
                for qubit_idx in range(len(z_coords)):
                    if z_coords[qubit_idx] == 0 and x_coords[qubit_idx] == 0:
                        operator.append("1")
                    elif z_coords[qubit_idx] == 0 and x_coords[qubit_idx] != 0:
                        operator.append("X" + ("" if x_coords[qubit_idx] == 1 else str(x_coords[qubit_idx])))
                    elif z_coords[qubit_idx] != 0 and x_coords[qubit_idx] == 0:
                        operator.append("Z" + ("" if z_coords[qubit_idx] == 1 else str(z_coords[qubit_idx])))
                    else:
                        operator.append("Z" + ("" if z_coords[qubit_idx] == 1 else str(z_coords[qubit_idx])) + "X" + ("" if x_coords[qubit_idx] == 1 else str(x_coords[qubit_idx])))
                row.append(operator)
            self.table.append(row)



    def print(self):
        for row in self.table:
            for operator in row:
                print(" ".join(operator) + "\t\t", end = "")
            print("\n")

