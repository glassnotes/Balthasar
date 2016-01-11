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

    def __init__(self, f, curve):
        # Set some obvious parameters
        self.p = f.p
        self.n = f.n
        self.dim = f.dim
        self.field = f # Keep a copy of the finite field to do math!
        self.curve = curve

        self.square = np.zeros((self.dim, self.dim), dtype=np.int)
        for i in range(self.dim):
            for j in range(self.dim):
                if self.n > 1:
                    self.square[i][j] = (self.field[j] + self.field.evaluate(self.curve, self.field[i])).prim_power
                else:
                    self.square[i][j] = (self.field[j] + self.field.evaluate(self.curve, self.field[i])).exp_coefs[0]


    def print(self):
        print(self.square)
