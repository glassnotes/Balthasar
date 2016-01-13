# Balthasar
A Python package for generating and manipulating mutually unbiased bases. Balthasar is currently under development and is as a result may not be fully functional. If you find any bugs, don't hesitate to contact me at odimatte@uwaterloo.ca.

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

### MUBs

To generate a set of MUBs in dimension _p_<sup>_n_</sup>, you must simply provide a finite field in this dimension:
```
from pynitefields import *
from balthasar import *

f = GaloisField(2, 3, [1, 1, 0, 1])  # Finite field of order 8
f.to_sdb([3, 5, 6])
dim8_mubs = MUBs(f)
```

It is best to use a self-dual, or almost self-dual basis to construct your MUB table. Using the
polynomial basis is not inherently a problem, but the program will produce a warning if you do so.


The created table of MUB operators can be viewed using
```
dim8_mubs.print()
```


### Curves
Balthasar contains a separate class Curve for manipulating curves over the finite field. 
A curve _b_(_a_) = _c_<sub>0</sub> + _c_<sub>1</sub> _a_ + ... + _c_<sub>_k_</sub> _a_<sup>k</sup>
where the _c_<sub>i</sub> are elements of the finite field is turned into a Curve object using a list
of the form [_c_<sub>0</sub>, _c_<sub>1</sub>, ..., _c_<sub>_k_</sub>] as in the following example:
```
f = GaloisField(2, 2, [1, 1, 1])
c = Curve([0, f[1], f[3], f[2]], f)
c.print() # Will print the curve in polynomial form
c.print(True) # Will print the points in the curve as tuples (repr as powers of the primitive element)
```
It is necessary to specify the field because we can have curves which simply contain constant coefficients:
```
f = GaloisField(5)
c = Curve([0, 1, 3], f) # [0, 1, 3] could be a curve in many different fields
```


### Latin squares

To generate a Latin square, one must simply pass it a curve over some finite field.
```
f = GaloisField(7)
c = Curve([0, f[1]], f)
l = LatinSquare(c)
l.print() # Prints the Latin square of order 7
```

### Advanced functionality
MUB operators can be viewed as living on rays in discrete phase space. Thus, one can provide a set 
of such rays to the constructor, as a list of objects of type curves.  You must provide a legitimate set of
d + 1 curves to do this - tread carefully.
If no set of rays is provided, the MUBs will be built from the set
of Desarguesian curves (linear curves). See the section on Curves below.
```
f = GaloisField(2, 2, [1, 1, 1])
rays = [Curve([0, el], f) for el in f]
rays.append(Curve([0, el], f, True))
some_mubs = MUBs(gf, rays)
```
=============================================

### Future features
- Explicit matrix forms of all operators using numpy
- Computation of mutually unbiased vectors from the operator table 
- Implementation of unitary transformations of the MUB table and associated curves
- Computation and transformations of associated Latin squares
- Coarse graining for incomplete tomography (research in progress!) 
- Wigner functions, point operators, and plotting
