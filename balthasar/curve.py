from pynitefields import *

class Curve():
    """ Class to hold all points in a curve.
    """
    def __init__(self, field, coefs, reverse=False):
        """
        Parameters:
        """
        self.field = field
        self.coefs = coefs
        self.form = "beta" # Default curve in form beta (alpha)
        self.is_ray = False # Does it pass through 0, 0
        self.points = [] # List of points as tuples

        if reverse == True: # Curve in form alpha (beta)
            self.form = "alpha"

        if type(coefs[0]) is int:
            if coefs[0] == 0:
                self.is_ray = True
        else:
            if coefs[0] == self.field[0]:
                self.is_ray = True

        for el in field:
            if self.form == "beta":
                self.points.append( (el, field.evaluate(coefs, el)) )
            elif self.form == "alpha":
                self.points.append( (field.evaluate(coefs, el), el) )


    def __iter__(self):                                                                                                                           
        """ Allow the user to iterate over the curve point by point """
        return iter(self.points)


    def print(self, as_points = False):
        if as_points == True: # Print the curve as a list of points
            for point in self.points:
                print(str(point[0].prim_power) + ", " + str(point[1].prim_power))
        else: # Print the curve as a polynomial
            for i in range(len(self.coefs)):
                if type(self.coefs[i]) == int: # Integers
                    if self.coefs[i] == 0: # Ignore 0 terms
                        continue
                    if self.coefs[i] == 1 and i == 0: # Only print 1 if it's the constant
                        print(str(self.coefs[i]), end="")
                else: # Field elements
                    if self.coefs[i] == self.field[0]: # Ignore 0 terms
                        continue
                    if (self.coefs[i].prim_power == 1):
                        print("s", end="")
                    else:
                        print("s" + str(self.coefs[i].prim_power), end="")

                if i > 0:
                    if i == 1:
                        print(" x", end="")
                    else:
                        print(" x" + str(i), end="")
                if i != len(self.coefs) - 1:
                    print(" + ", end = "")
        print("")

        
        
