Tutorial
************************************

Installation
====================================

Balthasar is written in Python 3. The following packages are required:

* PyniteFields_
* numpy
* matplotlib, libffi, cairocffi 

.. _PyniteFields: https://github.com/glassnotes/PyniteFields

Balthasar can be installed by running::

    python setup.py install

in the main directory of the program.

To generate a webpage version of this Sphinx documentation, one can run::

    make html

from the doc directory. Eventually this documentation will live somewhere 
online.

General functionality
====================================

Below we will go through some of the simple functionality of Balthasar:
making tables of MUB operators, Latin squares, and discrete Wigner Functions.
References for the underlying math can be found on the 
:ref:`references page <references>`.

Mutually unbiased bases
-------------------------------------
To generate a set of MUBs in dimension :math:`p^n`, you must simply provide 
a finite field expressed in the self-dual basis.::

    from pynitefields import *
    from balthasar import *

    # Construct a finite field and convert to self-dual basis (dimension 8 here)
    f = GaloisField(2, 3, [1, 1, 0, 1])  
    f.to_sdb([3, 5, 6])  

    # Construct the MUB by feeding it the field
    dim8_mubs = MUBs(f)

The created table of MUB operators can be viewed using print::

    dim8_mubs.print()

For large dimensions, constructing the MUBs can be quite intensive; by default,
matrix forms of all the operators are computed. To speed up initialization, we
can tell Balthasar not to compute the matrix, and just handle the operators
as names and strings. ::

    f256 = GaloisField(2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1])
    f256.to_sdb([5, 18, 30, 44, 106, 135, 147, 249])

    dim256_mubs = MUBs(f256) # Head over to reddit, this'll take a while
    dim256_mubs = MUBs(f256, matrix = False) # Runs in reasonable time

Unfortunately, MUBs constructed without matrix representations cannot be
used to numerically evaluate any Wigner functions. They can, however, be 
used to coarse-grain Wigner functions and compute new surviving displacement
operators. Balthasar will produce a warning when you construct MUBs without
matrices.

Note: currently, computation of the explicit mutually unbiased vectors is
not implemented. Apologies for the inconvenience; this is certainly on my 
list of things to do!

Discrete Wigner functions
-------------------------------------

One of the main purposes for writing Balthasar was to compute Wigner functions 
in discrete phase space. A Wigner function can be generated using a table of 
MUBs expressed in the self-dual basis. A Wigner function can be computed 
from either a state vector or a density matrix. Here we have an example in
dimension 4. ::

    # Initialize a field
    f = GaloisField(2, 2, [1, 1, 1])
    f.to_sdb([1, 2])

    # Construct MUB table
    dim4_mubs = MUBs(f)

    # Construct Wigner function framework by passing it MUBs
    my_wf = WignerFunction(dim4_mubs)

    # An arbitrary quantum state.
    state = (1.0 / math.sqrt(2)) * np.array([[1, 0, 0, 1]])

    # Compute the wf of state (returns a numpy array)
    wf_bell = my_wf.compute_wf(state)

    # Plot the Wigner function of state
    my_wf.plot(state)


Coarse-grained Wigner functions
--------------------------------------
Coarse-graining is the result of a foray into research on incomplete quantum
tomography. Our research is presented in [DMSSLG_]. 
A coarse-grained Wigner function is constructed from a fine-grained
one by aggregating the kernel (point) operators over a set of cosets of
the finite field. Much of the procedure for coarse graining has been automated,
though there are some tunable parameters. 

Let us begin with a simple set up: ::

    # We will coarse grain a dim 16 system to a dim 4 one
    f16 = GaloisField(2, 4, [1, 1, 0, 0, 1])
    f16.to_sdb([3, 7, 12, 13])

    f4 = GaloisField(2, 2, [1, 1, 1])
    f4.to_sdb([1, 2])

    m = MUBs(f16)

    fine_wf = WignerFunction(m)
    
There are two coarse-graining schemes we focus on. The first works in the 
general case; a basis is chosen for the big field with respect to the small 
field. This is the polynomial basis by default, but it is also possibly
to manually specify one using the optional 'basis' argument. ::

    # Coarse grain in general
    coarse_wf = CoarseWignerFunction(fine_wf, f4, mode='general')

For square dimensions, it is also possible to construct cosets using the 
copy of the small subfield within the big field. ::

    # Coarse grain using the subfield as the first coset
    coarse_wf = CoarseWignerFunction(fine_wf, f4, mode='subfield')

The class CoarseWignerFunction inherits from WignerFunction, so it is possible
to compute and plot Wigner functions like normal. ::

    from math import sqrt
    state = (1.0/sqrt(2))*np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])

    coarse_wf.compute_wf(state)
    coarse_wf.plot(state)

Finally, as per [DMSSLG_], coarse-graining also produces a subset of the displacement
operators which we purport can be used as measurement settings in an incomplete 
tomography scheme. These operators are accessible via the variable coarse_D in the CoarseWignerFunction ::

    print(coarse_wf.coarse_D)


Advanced functionality
=======================================

The remaining section of the tutorial pertains more to working with the 
structure of the underlying phase space than to the structures built on top
of it. 

In general, the set of 'standard' MUBs is associated to a bundle of lines, called rays,
of the form :math:`{\beta} = \lambda {\alpha}`. These are also called 
Desarguesian curves. However, there exist other sets of MUBs associated to 
different sets of curves. In most cases, these are unitarily equivalent to 
those of the Desarguesian set, but they may have different entanglement structures.

In what follows, we will discuss how to generate MUB tables using different
sets of curves. For these curves we can also plot their striations (sets of
parallel lines), and generate their associated Latin squares.

Curves
---------------------------------------
Balthasar contains a separate class Curve for manipulating curves over 
the finite field. A curve :math:`{\beta}({\alpha}) = c_0 + c_1 {\alpha} +
\cdots + c_k {\alpha}^k` where the :math:`c_i` are elements of the finite field 
is turned into a Curve object using a list of the form 
[:math:`c_0`, :math:`c_1`, ..., :math:`c_k`] as in the following example: ::

    f = GaloisField(2, 2, [1, 1, 1])
    c = Curve([0, f[1], f[3], f[2]], f)

    c.print() # Will print the curve in polynomial form
    c.print(True) # Will print the points in the curve as tuples 

Note that it is necessary to specify the field over which the curve is defined. 
This is because we can simplify cases where the coefficients are just integers: ::

    f = GaloisField(5)
    c = Curve([0, 1, 3], f) 

In general, curves are represented in phase space in the form :math:`{\beta} = 
f({\alpha})`. However, it is also possible to express a curve in the form
:math:`{\alpha}=f({\beta})` by passing an extra argument to function. ::

    f = GaloisField(5)
    cba = Curve([0, 1, 3], f)       # b = a + 3 a 
    cab = Curve([0, 1, 3], f, "alpha") # a = b + 3 b


Striations
-----------------------------------
Striations are the partitions of the affine plane into groups of parallel lines. 
They are used to build Latin squares and MUBs, and also to compute the point 
operators in discrete phase space for the Wigner function under the 
quantum net WF formulation found in [GHW_]. We are no longer using this formalism, but the 
striations are nevertheless useful to see, in particular when coarse-graining
Wigner Functions (where lines from the same striation are bundled
together and turned into thick lines). 

The set of striations of the rays :math:`{\beta} = f({\alpha})` can be generated using 
the code snippet below. The striations are stored as a list, with
the slopes in order of powers of the primitive element of the field, with the infinite slope last. 

You can view a striation graphically by using the plot function and passing 
in an index representing the desired power of the primitive element (and -1 for the vertical lines). ::

    f = GaloisField(2, 2, [1, 1, 1])
    s = Striations(f)
    s[0] # Rays
    s.plot()  # Graphically see the rays
    s.plot(2) # Plot the striation with slope x^2, x is the primitive element of f


Latin squares
--------------------------------------
Latin squares can be contructed from non-degenerate curves over finite fields.
'Non-degenerate' means that the curve is something called a permutation
polynomial, i.e. putting the entire field through the curve gives us the field
back in a permuted order.

As MUBs can be associated with sets of non-degenerate curves (which are also
additive and commutative), we can consider that some MUBs can be associated
with complete sets of Latin squares. These Latin squares have a special property,
that of being a complete, mutually orthogonal set. Some unitary transformations 
on these MUBs sometimes lead to a new set of mutually orthogonal Latin squares
which is isomorphic to the first. These relationships are discussed in detail
in previous work, [GDMKdG_].

To generate a Latin square in Balthasar, one must simply pass it a curve 
over some finite field. ::

    f = GaloisField(7)
    c = Curve([0, f[1]], f)
    l = LatinSquare(c)
    l.print() # Prints the Latin square of order 7



MUBs and curves
--------------------------------------
By default, MUBs will be constructed with the set of Desarguesian curves.
However, we can specify a set of d + 1 curves with which to produce MUBs.
We show here an example in dimension 4. The set of curves is taken from
[KRBSS07_]. ::

    f = GaloisField(2, 2, [1, 1, 1])

    c1 = Curve([0, 0, f[3]], f)    # beta = alpha^2
    c2 = Curve([0, f[3], f[3]], f) # beta = alpha + alpha^2
    c3 = Curve([0, f[1], f[3]], f) # beta = sigma alpha + alpha^2
    c4 = Curve([0, f[2], f[3]], f) # beta = sigma^2 alpha + alpha^2
    c5 = Curve([0, 0], f, "alpha")    # alpha = 0 
    
    curves = [c1, c2, c3, c4, c5]

    some_mubs = MUBs(f, curves)


.. _KSSdG: http://iopscience.iop.org/article/10.1088/0305-4470/38/12/015/meta
.. _RBKSS: http://journals.aps.org/pra/abstract/10.1103/PhysRevA.72.062310
.. _KRBSS09: http://www.sciencedirect.com/science/article/pii/S0003491608001541
.. _KRBSS07: http://iopscience.iop.org/article/10.1088/1751-8113/40/14/014/meta
.. _GHW: http://journals.aps.org/pra/abstract/10.1103/PhysRevA.70.062101
.. _GDMKdG: http://iopscience.iop.org/article/10.1088/1751-8113/47/43/435303/meta
.. _DMSSLG: 
