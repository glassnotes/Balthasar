from pynitefields import *
from balthasar.curve import Curve

class Striations():
    """ A class to keep the set of linear striations in a given dimension
        for use in the Wigner function and doing various other things.
        It'll just be cleaner to have them in a separate class, and might
        be useful for other things later too.

        Striations are essentially stored as a set of Curves over the provided
        finite field.
    """
    def __init__(self, field):
        self.field = field
        self.striations = []

        # Handle the horizontal striation first
        horizontal_striation = []
        for intercept in field:
            curve = Curve([intercept, self.field[0]], self.field, True)
            horizontal_striation.append(curve)
        self.striations.append(horizontal_striation)

        # Take care of the rest of the striations 
        for slope in field:
            new_striation = []
            for intercept in field:
                    curve = Curve([intercept, slope], self.field)
                    new_striation.append(curve)
            self.striations.append(new_striation)
    
        # Store the rays as a set for safekeeping
        rays = [s[0] for s in self.striations]


    def print(self, as_points = False):
        print("============================")
        for s in self.striations:
            print("Ray: ")
            s[0].print()

            for curve in s:
                curve.print(as_points)
                print("")
            print("============================")
        


            
