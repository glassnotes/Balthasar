from pynitefields import *
from balthasar.striations import Striations
import numpy as np

class WignerFunction():
    """ Class to store and plot a discrete Wigner function.

        Parameters:
        ===========
        field - The finite field over which the Wigner function is defined
        mubs - The MUBs associated with this Wigner function
        quantum_net - Association of striations to density operators
        striations - The Striation
        point_operators - Operators at each point in discrete phase space
        matrix - The matrix which contains the probability values for each Wigner function
    """

    def __init__(self, mubs):
        """ Initialize the Wigner function and compute the point operators. """
        self.field = mubs.field
        self.mubs = mubs
        self.striations = Striations(self.field)
        self.quantum_net = self.compute_quantum_net()
        self.point_operators = self.compute_point_operators(self.striations, self.quantum_net)    

         
    def compute_quantum_net(self):
        """ Compute the quantum net of the Wigner function, i.e. the association of
            lines with quantum states (as density operators).

            We associate every ray with the all +1 combination of the generators of
            the row of MUB operators associated to that curve. 
        """
        # Get the table of generating operators and their matrix representation
        gen_ops, gen_mats = self.mubs.compute_generators() 

        # Horizontal shifts correspond to the Z-only operators (b = 0, so only alpha coordinate varies)
        # Vertical shifts correspond to the X-only operators (a = 0, so only beta coordinate varies)
        trans_hor, trans_vert = self.mubs.matrix_table[0], self.mubs.matrix_table[-1]
        
        quantum_net = []    

        # Loop through every striation
        for str_idx in range(len(gen_ops)):
            next_striation_net = [] 
            for line_idx in range(len(self.striations[0])): # For every line in the striation
                if line_idx == 0: # If it's the ray... 
                    gen_sum = (1.0 / self.field.dim) * np.sum(gen_mats[str_idx], 0) # Sum the generators
                    next_striation_net.append(gen_sum)
                else: # Otherwise, transform according to one of the translation operators
                    transformation = np.zeros((self.field.dim, self.field.dim))
                    if str_idx == (len(gen_ops) - 1): # Vertical striation must be handled separately and shifted using Z operators
                        transformation = trans_hor[line_idx - 1] # -1 accounts for the fact that the identity is technically the 0th trans
                    else: # All other striations shifted vertically horizontally using the X operators
                        transformation = trans_vert[line_idx - 1]
                    transformed_operator = np.dot(np.dot(transformation, next_striation_net[0]), np.asmatrix(transformation).H)
                    next_striation_net.append(transformed_operator)
            quantum_net.append(next_striation_net)
                        
        return quantum_net        


    def compute_point_operators(self, striations, net):
        """ Use the quantum net and the striations to compute the point operators
            for the Wigner function. 
        """
        # The Wigner function is going to be "upside-down" for now
        # A point operator is the sum of line operators through that point minus the identity
        point_operators = []

        for beta in self.field: # Fix the vertical axis
            point_ops_this_row = []
            for alpha in self.field: # Move across horizontally
                point = (alpha, beta)
                point_op = np.zeros((striations.field.dim, striations.field.dim))
                # Every point will appear once in each striation
                for str_idx in range(self.field.dim + 1):
                    for curve_idx in range(self.field.dim):
                        if point in striations[str_idx][curve_idx]: # Found the point
                            point_op = point_op + net[str_idx][curve_idx] # Add the operator
                            continue # Don't do more work than we have to, this is already hideous
                point_op = point_op - np.eye(striations.field.dim) # Subtract the identity
                point_ops_this_row.append(point_op)
            point_operators.append(point_ops_this_row)
        return point_operators
                            
                                

    def compute_wf(self, state):
        """ Compute the probabilities in the Wigner function for a given state.
            Input: state, a numpy array representing either a ket or a density matrix.
        """
        W = np.zeros(( self.field.dim, self.field.dim )) # Holds result

        # Don't discriminate - allow the user to submit either a ket state vector
        # or a density operator. If the state is a ket, just turn it into density operator.
        if state.shape[0] == 1:
            state = np.outer(state, np.conj(state))

        for beta in range(self.field.dim):
            for alpha in range(self.field.dim):
                W[beta][alpha] = 1.0 / self.field.dim * np.trace(np.dot(state, self.point_operators[beta][alpha]))

        return W


    def plot(self, state):
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
            
