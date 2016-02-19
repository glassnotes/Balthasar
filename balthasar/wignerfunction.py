from pynitefields import *
import numpy as np

class WignerFunction():
    """ Class to store and plot a discrete Wigner function.

        Parameters:
        ===========
        field - The finite field over which the Wigner function is defined
        point_operators - Operators at each point in discrete phase space
        matrix - The matrix which contains the probability values for each Wigner function
        is_coarse - Optional parameter, used for specifying whether the function is coarse grained
                    (true) or fine grained (false, i.e. a normal Wigner function)
    """

    def __init__(self, field):
        """ Initialize the Wigner function and compute the point operators. """


    def compute_wf(state):
        """ Compute the probabilities in the Wigner function for a given state.
            Input: state, a numpy array representing either a ket or a density matrix.
        """
        dim = len(point_ops)
        W = np.zeros(( int(np.sqrt(dim)), int(np.sqrt(dim)) )) # Holds result

        # Don't discriminate - allow the user to submit either a ket state vector
        # or a density operator. If the state is a ket, just turn it into density operator.
        if state.shape[0] == 1:
            state = np.outer(state, np.conj(state))

        for point in point_ops:
            W[point[0]][point[1]] = 1.0 / dim * np.trace(np.dot(state, point_ops[point]))

        return W


    def plot(state):
        """
        Plot the Wigner function of a given state.
        """
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        W = compute_wf(state)

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
        plt.gca().set_zlim([0, 0.25])
            
