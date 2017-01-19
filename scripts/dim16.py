""" A sample script to try out coarse-graining in dimension 16 
    down to dimension 4. Feel free to run, modify, and take
    snippets of it to use for other things!

    Author: Olivia Di Matteo, 2017
"""

from pynitefields import *
from balthasar import *

import numpy as np
from math import sqrt

# Generate the finite fields F16 and F4
f16 = GaloisField(2, 4, [1, 1, 0, 0, 1])
f16.to_sdb([3, 7, 12, 13])
f4 = GaloisField(2, 2, [1, 1, 1])
f4.to_sdb([1, 2])

# Generate the MUBs
m = MUBs(f16)

# Generate the Wigner function
wf = WignerFunction(m)

# Create a state as a numpy vector 
# In this case we choose two Bell states tensored together
# Make sure it's normalized!
s = np.zeros((1, 16))
s[0][0] = 1.0 /(2)
s[0][3] = 1.0/(2)
s[0][12] = 1.0/(2)
s[0][-1] = 1.0 / (2)

# Plot the Wigner function in a new window
wf.plot(s)

# Or plot it in a file
# wf.plot(s, 'wf_doublebell.png')

# Create a coarse-grained Wigner function using the 'general' method 
cwf = CoarseWignerFunction(wf, f4) 

# Print out the matrix verion and plot it
print(cwf.compute_wf(s))
cwf.plot(s) 

# Print out the surviving displacement operators for incomplete tomo
from pprint import pprint
pprint(cwf.coarse_D)

# Coarse grain again by cosetting w.r.t. the subfield copy of F4 in F16 
cwf_sub = CoarseWignerFunction(wf, f4, mode = 'subfield')
cwf_sub.plot(s)
pprint(cwf_sub.coarse_D)
