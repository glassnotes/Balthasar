""" A sample script to try out coarse-graining in dimension 256 
    down to dimension 16. Feel free to run, modify, and take
    snippets of it to use for other things!

    As dimension 256 is rather large, we show here an example
    where we do not compute the full MUB matrix table or 
    values of any Wigner functions, but rather just coarse-grain
    to obtain a subset of the displacement operators.

    The obtain of this scripts is a series of lists indexed by 
    a number - these numbers are the slopes of the rays in which
    these surviving MUBs should be put in the coarse-grained phase-space.
    Note that there are horizontal/vertical cases, as well as 
    precisely 15 intermediate cases (those corresponding to the
    elements of the subfield F16 in F256). See our recent
    paper on the topic for more details.

    Author: Olivia Di Matteo, 2017
"""

from pynitefields import *
from balthasar import *

import numpy as np
from math import sqrt

# Generate the finite fields F16 and F4
f256 = GaloisField(2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1])
f256.to_sdb([5, 18, 30, 44, 106, 135, 147, 249])
f16 = GaloisField(2, 4, [1, 1, 0, 0, 1])
f16.to_sdb([3, 7, 12, 13])

# Generate the MUBs
# Will print out a warning because there are no matrices.
m = MUBs(f256, matrix = False)

# Generate the Wigner function
wf = WignerFunction(m)

# Create a coarse-grained Wigner function using the 'general' method 
cwf = CoarseWignerFunction(wf, f16) 

# Print out the surviving displacement operators for incomplete tomo
from pprint import pprint
pprint(cwf.coarse_D)

# Coarse grain again by cosetting w.r.t. the subfield copy of F16 in F256
cwf_sub = CoarseWignerFunction(wf, f16, mode = 'subfield')
pprint(cwf_sub.coarse_D)
