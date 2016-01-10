# Balthasar
A Python package for generating and manipulating mutually unbiased bases

Dependencies

===================================================

To run the routines in Balthasar, you will need to have installed my PyniteFields package
(available at github.com/glassnotes/PyniteFields).

You will also need numpy.

===================================================

To generate a set of MUBs in dimension p^n is done by providing a finite field in this dimension:
```
from pynitefields import *
from balthasar import *

gf = GaloisField(2, 3, [1, 1, 0, 1])  # Finite field of order 8
gf.to_sdb([3, 5, 6])
dim8_mubs = MUBs(gf)
```

It is best to use a self-dual, or almost self-dual basis to construct your MUB table. Using the
polynomial basis is not inherently a problem, but the program will produce a warning if you do so.


The created table of MUB operators can be viewed using
```
dim8_mubs.print()
```

MUB operators can be viewed as living on rays in discrete phase space. Thus, one can provide a set 
of such rays to the constructor. If no set of rays is provided, the MUBs will be built from the set
of Desarguesian curves (linear curves) of the form $$ \beta = \lambda \alpha $$.
```
some_mubs = MUBs(gf, rays)
```

