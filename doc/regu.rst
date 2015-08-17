.. _regu:

Test Matrices for Regularization Method
----------------------------------------

A Fredholm integral equation of the first kind (in 1-dimensional) can
be written as
 
.. math::

   \int_{0}^1 K(s,t) f(t) dt = g(s), \quad 0 \leq s \leq 1,

where :math:`g` and :math:`K` are known functions and :math:`f` is the
unknown solution. This is a classical example of a linear ill-posed
problem, i.e., an arbitrary small perturbation of the data can cause
an arbitrarily large perturbation of the solution. For example, in
computerized tomography, :math:`K` is an X-ray source, :math:`f` is
the object being scanned, and :math:`g` is the measured damping of the
X-rays. The goal here is to reconstruct the scanned object from
information about the locations of the X-ray sources and measurements
of their damping. 

After discretizations (by the quadrature method or the Galerkin
method), we obtain a linear system of equations :math:`Ax=b`. The
generated test problems in Matrix Depot has type ``RegProb``, which is
defined as::

  immutable RegProb{T}
    A::Matrix{T}  # matrix of interest
    b::Vector{T}  # right-hand side
    x::Vector{T}  # the solution to Ax = b
  end

