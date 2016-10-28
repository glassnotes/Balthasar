from pynitefields import *
from balthasar.striations import Striations
from functools import reduce
import numpy as np
import pprint as pp

class WignerFunction():
    """ Class to store and plot a discrete Wigner function.

        Parameters:
        ===========
        field - The finite field over which the Wigner function is defined
        dim - The dimension of the system 
        mubs - The MUBs associated with this Wigner function
        D - The dictionary of displacement operators
        kernel - Operators at each point in discrete phase space ('point ops')
    """

    def __init__(self, mubs):
        """ Initialize the Wigner function and compute the point operators. """
        self.field = mubs.field
        self.dim = self.field.dim
        self.mubs = mubs
        self.D = mubs.D
        self.kernel = self.compute_kernel()    

         
    def compute_kernel(self):
        """ Compute the 'kernel' of the Wigner function, i.e. the set of 
            point operators.
        """
        kernel = {} 

        # Computation of the kernel can be accomplished by computing the
        # value at point (0, 0), then translating it using the D operators,
        # i.e.        w(a, b) = D(a, b) w(0, 0) D(a, b)^\dag
          
        kernel_00 = np.zeros((self.dim, self.dim), dtype=np.complex_)
        if self.field.p != 2:
            for key in self.D.keys():
                kernel_00 = kernel_00 + (self.D[key][0].eval() * self.D[key][1])
            kernel_00 = kernel_00 / self.dim
            kernel[(self.field[0], self.field[0])] = kernel_00
        else:
            for key in self.D.keys():
                kernel_00 = kernel_00 + (self.D[key][0] * self.D[key][1])
            kernel_00 = kernel_00 / self.dim
            kernel[(self.field[0], self.field[0])] = kernel_00

        # Compute the rest of the points by translating the first one
        for a in self.field:
            for b in self.field:
                if a == self.field[0] and b == self.field[0]:
                    continue # Don't set the 0 case again
                # Compute the displacement operator with phase included
                if self.field.p == 2: # For p = 2 this is a number
                    dab = self.D[(a, b)][0] * self.D[(a, b)][1]
                else: # For p != 2 we need to evaluate the power of pth root 
                    dab = self.D[(a, b)][0].eval() * self.D[(a, b)][1]

                dab_dag = np.asmatrix(dab).getH()

                kernel[(a, b)] = np.dot(dab, np.dot(kernel_00, dab_dag))

        return kernel

                            
                                
    def compute_wf(self, state):
        """ Compute the probabilities in the Wigner function for a given state.
            Input: state, a numpy array representing either a ket or a 
            density matrix.
        """
        W = np.zeros((self.dim, self.dim)) # Holds result

        # Don't discriminate - allow the user to submit either a ket vector
        # or a density operator. If it's a ket, switch to a density operator.
        if state.shape[0] == 1:
            state = np.outer(state, np.conj(state))

        sorted_els = sorted(self.field.elements)

        for a in self.field:
            for b in self.field:
                mat = np.trace(np.dot(state, self.kernel[(a, b)]))
                a_coord = sorted_els.index(a)
                b_coord = sorted_els.index(b)
                #W[a.prim_power][b.prim_power] = (1.0 / self.dim) * mat
                W[a_coord][b_coord] = (1.0 / self.dim) * mat
                

        return W

    
    def plot_mat(self, state):
        import matplotlib.pyplot as plt
        W = self.compute_wf(state)
        plt.matshow(W)
        plt.show()


    
    def plot(self, state):
        """
        Plot the Wigner function of a given state.
        """
        
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        W = self.compute_wf(state)

        # This code was written in the summer and frankly I still don't 
        # fully understand, or at this point remember, why it works >.<
        fig = plt.figure()

        ax1 = fig.add_subplot(211, projection='3d')

        xpos = []
        ypos = []
        zpos = []

        for x in range(0, len(W)):
            for y in range(0, len(W)):
                xpos.append(x)
                ypos.append(y)
                zpos.append(0)

        dx = np.ones(len(xpos))
        dy = np.ones(len(xpos))
        dz = W.flatten()

        max_colour = dz.max()
        min_colour = dz.min()
        colours = []

        for point in dz:
            if point < 0:
                colours.append(( 1 - (point / min_colour), 1 - (point / min_colour), 1))
            elif point == 0:
                colours.append((1, 1, 1))
            else:
                colours.append((1, 1 - (point / max_colour), 1 - (point / max_colour)))

        ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colours)

        plt.gca().invert_xaxis()
        plt.xticks(range(0, len(W)))
        plt.yticks(range(0, len(W)))
        plt.gca().set_zlim([0, dz.max()])
        plt.show()

            
