# Balthasar
A Python package for generating and manipulating mutually unbiased bases.

## Dependencies

- [PyniteFields](https://github.com/glassnotes/PyniteFields)
- numpy

===================================================

## Installation
In the main directory of the program, you can install Balthasar into your version of Python 3 using the following command:

```
python3 setup.py install
```

## Usage

To generate a set of MUBs in dimension _p_<sup>_n_</sup>, you must simply provide a finite field in this dimension:
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
of Desarguesian curves (linear curves). 
```
some_mubs = MUBs(gf, rays)
```

