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
        self.dim = field.dim
        self.striations = []

        # Striations of the form beta = lambda alpha + gamma
        # Do these first so that the slope lambda corresponds to the list element 
        for slope in field:
            new_striation = []
            for intercept in field:
                    curve = Curve([intercept, slope], self.field)
                    new_striation.append(curve)
            self.striations.append(new_striation)

        # Handle the vertical striation (infinite slope) last
        # This way it is easily accessible as striation -1
        horizontal_striation = []
        for intercept in field:
            curve = Curve([intercept, self.field[0]], self.field, True)
            horizontal_striation.append(curve)
        self.striations.append(horizontal_striation)

        # Store the rays as a set for safekeeping
        rays = [s[0] for s in self.striations]

        

    def __iter__(self):
        """ Allow the user to iterate through all the striations. """
        return iter(self.striations)


    def __getitem__(self, index):
        """ Grab a striation. Striation -1 corresponds to the horizontal striation,
            all other striations indicated by there slope as a power of the 
            primitive field element.
        """
        if (index >= -1) and (index <= self.dim + 1):
            return self.striations[index]
        else:
            print("Error, element out of bounds.")


    def get_rays(self):
        """ Get the set of rays only. """
        return [s[0] for s in self.striations]


    def print(self, as_points = False):
        """ Print out all the striations either as equations or as sets of points. """
        print("============================")
        for s in self.striations:
            print("Ray: ")
            s[0].print()

            for curve in s:
                curve.print(as_points)
                print("")
            print("============================")
