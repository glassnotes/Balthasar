from pynitefields import *
import numpy as np

class LatinSquare():
    """ Class to hold a Latin square.
        
        Member variables:
        p - Prime dimension
        n - Power of prime 
        dim - The dimension of the system p^n
        field - The finite field used to build this square 
        curve - An object of type Curve stored in the field
        square - The actual square (stored as a numpy array)
    """

    def __init__(self, curve):
        # Set some obvious parameters
        self.field = curve.field # Keep a copy of the finite field to do math!
        self.dim = self.field.dim
        self.curve = curve

        self.square = np.zeros((self.dim, self.dim), dtype=np.int)

        for row in range(self.dim):
            for column in range(self.dim):
                row_el = curve[row][1]
                col_el = curve[column][0]
                self.square[row][column] = (row_el + col_el).prim_power


    def print(self):
        print(self.square)
