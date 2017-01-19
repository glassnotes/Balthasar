""" A sample script to try out coarse-graining in dimension 9 
    down to dimension 3. Feel free to run, modify, and take
    snippets of it to use for other things!

    Qudit coarse-graining currently works with a few hiccups
    in terms of the plotting. I'll work on fixing these
    eventually, but the mechanics of coarse-graining itself work!

    Author: Olivia Di Matteo, 2017
"""

from pynitefields import *
from balthasar import *

import numpy as np
from math import sqrt

# Generate the finite fields F9 and F3
f9 = GaloisField(3, 2, [2, 1, 1])
# The self-dual basis in dimension 9 is *almost* self-dual 
f9.to_sdb([2, 4])

f3 = GaloisField(3)

# Generate the MUBs
m = MUBs(f9)

# Generate the Wigner function
wf = WignerFunction(m)

# Create a state as a numpy vector 
# Make sure it's normalized!
s = np.zeros((1, 9))
s[0][0] = 1.0 / sqrt(3)
s[0][4] = 1.0 / sqrt(3)

# Plot the Wigner function in a new window
wf.plot(s)

# Create a coarse-grained Wigner function using the 'general' method 
cwf = CoarseWignerFunction(wf, f3) 

# Print out the matrix verion and plot it
print(cwf.compute_wf(s))
cwf.plot(s) 

# Print out the surviving displacement operators for incomplete tomo
from pprint import pprint
pprint(cwf.coarse_D)

# Coarse grain again by cosetting w.r.t. the subfield copy of F3 in F9 
cwf_sub = CoarseWignerFunction(wf, f3, mode = 'subfield')
cwf_sub.plot(s)
pprint(cwf_sub.coarse_D)
