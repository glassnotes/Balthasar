from pynitefields import *
import numpy as np

class LatinSquare():
    """ Class to hold a Latin square.
        
        Member variables:
        p - Prime dimension
        n - Power of prime 
        dim - The dimension of the system p^n
        field - The finite field used to build this square 
        curve - Coefficients of the curve in the field 
        square - The actual square (stored as a numpy array)
    """

    def __init__(self, curve, f):
        # Set some obvious parameters
        self.p = f.p
        self.n = f.n
        self.dim = f.dim
        self.field = f # Keep a copy of the finite field to do math!
        self.curve = curve

        self.square = np.zeros((self.dim, self.dim), dtype=np.int)

        for row in range(self.dim):
            for column in range(self.dim):
                row_el = curve[row][1]
                col_el = curve[column][0]

                self.square[row][column] = (row_el + col_el).prim_power


    def print(self):
        print(self.square)
