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
            if coefs[0] == f[0]:
                self.is_ray = True

        for el in field:
            if self.form == "beta":
                self.points.append( (el, field.evaluate(coefs, el)) )
            elif self.form == "alpha":
                self.points.append( (field.evaluate(coefs, el), el) )

    def print(self):
        for point in self.points:
            print(str(point[0].prim_power) + ", " + str(point[1].prim_power))
        
        
