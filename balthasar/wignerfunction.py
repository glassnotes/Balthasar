#-*- coding: utf-8 -*-

from pynitefields import *

import numpy as np

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

        if self.mubs.matrices == False:
            print("Warning: Wigner Function constructed without matrices.")
            print("This Wigner Function cannot be used for plotting or \
                any numerical purposes, aside for computing coarse-grained \
                displacement operators.")
            return

        self.D = mubs.D
        self.kernel = self.compute_kernel()    

         
    def compute_kernel(self):
        """ Compute the 'kernel' of the Wigner function, i.e. the set of 
            what Wootters calls point operators.
        """
        kernel = {} 

        if self.mubs.matrices == False:
            print("Error, no matrices, cannot compute Wigner Function kernel.")
            return

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
        if self.mubs.matrices == False:
            print("Error, no matrices, cannot compute Wigner Function.")
            return

        W = np.zeros((self.dim, self.dim)) # Holds result

        # Don't discriminate - allow the user to submit either a ket vector
        # or a density operator. If it's a ket, switch to a density operator.
        if state.shape[0] == 1:
            state = np.outer(state, np.conj(state))

        # We sort the elements of the finite field so that Wigner function
        # Gets plotted in the "correct" order, i.e. the computational basis
        # states go 000, 001, 010, etc. and same for the +/- basis. This in 
        # a way imposes an 'order' on the field elements.
        sorted_els = sorted(self.field.elements)

        for a in self.field:
            for b in self.field:
                mat = np.trace(np.dot(state, self.kernel[(a, b)]))
                a_coord = sorted_els.index(a)
                b_coord = sorted_els.index(b)
                W[a_coord][b_coord] = (1.0 / self.dim) * mat

        return W

    
    def plot_mat(self, state):
        """ A simple matrix plot of the Wigner function. """
        import matplotlib.pyplot as plt
        W = self.compute_wf(state)
        plt.matshow(W)
        plt.show()


    
    def plot(self, state, filename=""):
        """ Compute and plot the Wigner function of a given state. 
            The parameter state can be either a vector or full 
            density matrix.
            
            If a filename is passed as a parameter, the plot function
            will output the plot to an eps file rather than display it.
        """
    
        if self.mubs.matrices == False:
            print("Error, no matrices, cannot plot Wigner function.")
            return

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

        # For standard Wigner functions, set the ticks to be basis states.
        # When the Wigner function gets created it is ordered in this way.
        if type(self) == WignerFunction:
            fmt_str = "#0" + str(self.field.n + 2) + "b"
            comp_basis = [format(x, fmt_str)[2:] for x in range(self.field.dim)]
            pm_basis = [x.replace('0', '+').replace('1', '\u2013') for x in comp_basis]
            ax1.set_xticklabels(pm_basis)
            ax1.set_yticklabels(comp_basis)
        else:
            plt.xticks(range(0, len(W)))
            plt.yticks(range(0, len(W)))
        
        plt.gca().set_zlim([0, dz.max() + 0.0001])

        plt.tight_layout()

        if filename == "":
            plt.show()
        else:
            plt.savefig(filename, dpi=1200, bbox_inches='tight')

