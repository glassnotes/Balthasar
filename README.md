# Balthasar
A Python 3 package for generating and manipulating mutually unbiased bases and 
various other related objects. Balthasar is under continuous development and as 
a result may not be fully functional. If you find any bugs, or have any questions,
don't hesitate to contact me at odimatte@uwaterloo.ca.

## Dependencies

- [PyniteFields](https://github.com/glassnotes/PyniteFields)
- numpy
- matplotlib, libffi, cairocffi (for plotting Wigner functions)

===================================================

## Installation
In the main directory of the program, you can install Balthasar into your version of Python 3 using the following command:

```
python3 setup.py install
```

## Usage

Below is a brief tutorial to get things started. Balthasar also comes with 
Sphinx documentation - this can be generated by
entering the doc directory and running, for example, 'make html'. I've included
a pdf in the main directory, and there is a more detailed tutorial in there.

### MUBs

To generate a set of MUBs in dimension _p_<sup>_n_</sup>, you must simply 
provide a finite field or order _p_<sup>_n_</sup> expressed in a self-dual basis. Currently Balthasar 
fully supports multi-qubit systems, with work on multi-qudit systems in progress.

```
from pynitefields import *
from balthasar import *

f = GaloisField(2, 3, [1, 1, 0, 1])  # Finite field of order 8
f.to_sdb([3, 5, 6])  # Convert to self-dual basis
dim8_mubs = MUBs(f)
```

Some possibilities for fields and bases are provided in the README for 
[PyniteFields](https://github.com/glassnotes/PyniteFields).

The created table of MUB operators can be viewed using

```
dim8_mubs.print()
```

### Wigner functions

One of the main purposes for writing Balthasar was to compute Wigner functions 
in discrete phase space. A Wigner function can be generated using a table of 
MUBs expressed in a self-dual basis. A Wigner function can be computed from 
either a state vector or a density matrix.

```
f = GaloisField(2, 2, [1, 1, 1])
f.to_sdb([1, 2])

dim4_mubs = MUBs(f)

my_wf = WignerFunction(dim4_mubs)

state = (1.0 / math.sqrt(2)) * np.array([[1, 0, 0, 1]])
wf_bell = my_wf.compute_wf(state) # Compute the WF of a state

my_wf.plot(state) # Plot the WF in 3d
```

### Coarse-graining the Wigner function

One of the primary reasons for putting together Balthasar was to
make it easy for researchers to make use of the procedure we present
in our recent paper [citation coming soon!]. A coarse Wigner function
can be created by specifying an existing fine-grained function, and a 
subfield. Coarse-grained functions and plots can then be computed in the
same manner as normal Wigner functions. Furthermore, a by-product of 
coarse-graining is a subset of the displacement operators; these are
stored and can be easily displayed as in the example below.

```
# Create a 'big' field, and a 'small', coarse field
field = GaloisField(2, 4, [1, 1, 0, 0, 1]) # Dimension 16
field.to_sdb([3, 7, 12, 13])

subfield = GaloisField(2, 2, [1, 1, 1]) # Dimension 4
subfield.to_sdb([1, 2])

mubs = MUBs(field)
wf = WignerFunction(mubs)

coarse_wf = CoarseWignerFunction(wf, subfield)

state = # Some 4-qubit state

coarse_wf_matrix = coarse_wf.compute_wf(state) # Compute the coarse WF for state
coarse_wf.plot(state) # Plot the coarse WF of state

print(coarse_wf.coarse_D) # Output the surviving coarse displacement operators

```


### Curves
Balthasar contains a separate class Curve for manipulating curves over the 
finite field. A curve 
_b_(_a_) = _c_<sub>0</sub> + _c_<sub>1</sub> _a_ + ... + _c_<sub>_k_</sub> _a_<sup>k</sup> 
where the _c_<sub>i</sub> are elements of 
the finite field is turned into a Curve object using a list of the form 
[_c_<sub>0</sub>, _c_<sub>1</sub>, ..., _c_<sub>_k_</sub>] as in the 
following example:

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

Finally, curves are generally thought of as living in discrete phase space, a two-dimensional grid _p_<sup>_n_</sup> x _p_<sup>_n_</sup> grid. In this space,
_a_ is the coordinate on the horizontal axis, and _b_ on the vertical. Usually we write our curves as _b_(_a_), however some curves must be specified 'backwards' in phase space, i.e. in the form _a_(_b_). This can be
achieved by passing the optional mode parameter "alpha" (and equivalently
"beta" for completeness).

```
c_b0 = Curve([0, 0], f) # Curve beta = 0, the horizontal line
c_a0 = Curve([0, 0], f, "alpha") # Curve alpha = 0, the vertical line
```

### Striations
Striations are the partitions of the affine plane into groups of parallel lines. 
They are used to build Latin squares (and can make for some very pretty 
graphics). The _rays_ (the straight lines passing through the origin) are an 
important set of curves, and are often used to construct MUBs, as they produce
a valid, complete set of MUBs in every dimension where such a set exists.

The set of striations can be generated using the code snippet below; the first
element of the striation list is the aforementioned rays.

You can view a striation graphically by using the plot function and passing 
in an index. Striations are indexed by their slope, from 0 through each field 
element to 1, and then the last (-1) striation is the vertical lines with 
infinite slope.

```
f = GaloisField(2, 2, [1, 1, 1])
s = Striations(f)
s[0] # Rays
s.plot()  # Graphically see the rays
s.plot(2) # Plot the striation with slope x^2, x is the primitive element of f
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
MUB operators can be constructed from curves in discrete phase space. By default,
Balthasar will construct MUBs based on rays, however other sets of MUB exist which
are defined in terms of a bundle of additive, commutative curves.
One can provide a set of such curves to the constructor, as a list of _d_ + 1
objects of type Curve.  These curves must satisfy the 
conditions listed in section 3.1 of 
[this paper](http://www.sciencedirect.com/science/article/pii/S0003491608001541). Currently there is no mechanism in place to check the validity of the curves; be careful! 

Below is an example which manually constructs the rays and passed them to the
MUB constructor.

```
f = GaloisField(2, 2, [1, 1, 1])
rays = [Curve([0, el], f) for el in f]
rays.append(Curve([0, el], f, "alpha"))
some_mubs = MUBs(gf, rays)
```
=============================================

### Future features
- Computation of mutually unbiased vectors from the operator table (this is 
  harder than it sounds)
- Implementation of unitary transformations of the MUB table and associated curves
- Transformations of Latin squares
- Other fun objects (e.g. unitary error bases, but also open to suggestions)
