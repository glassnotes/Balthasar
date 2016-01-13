from pynitefields import *
from balthasar.curve import Curve 

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
        if curves == []: # Default to Desarguesian curves
            self.curves.append(Curve(self.field, [0, f[0]], True)) # Horizontal curve alpha = 0
            for el in self.field: # Rest of the curves
                self.curves.append(Curve(self.field, [0, el]))
        else: # Curves specified by user
            if self.verify_curves(curves) == True:
                self.curves = curves 

        # Build the operator table
        self.table = []
    
        for curve in self.curves:
            curve.print()
            row = []

            for point in curve:
                op= []
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

            self.table.append(row)


    def verify_curves(self, curves):
        # TODO check that the properties of the provided curves are valid
        if len(curves) != self.dim + 1:
            print("Error, not enough curves provided.")
            return False
        return True


    def print(self):
        for row in self.table:
            for operator in row:
                print(" ".join(operator) + "\t\t", end = "")
            print("\n")

