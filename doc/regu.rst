.. _regu:

Test Problems for Regularization Method
=======================================

A Fredholm integral equation of the first kind (in 1-dimensional) can
be written as
 
.. math::

   \int_{0}^1 K(s,t) f(t) dt = g(s), \quad 0 \leq s \leq 1,

where :math:`g` and :math:`K` (called kernel) are known functions
and :math:`f` is the unknown solution. This is a classical example of
a linear ill-posed problem, i.e., an arbitrary small perturbation of
the data can cause an arbitrarily large perturbation of the
solution. For example, in computerized tomography, :math:`K` is an
X-ray source, :math:`f` is the object being scanned, and :math:`g` is
the measured damping of the X-rays. The goal here is to reconstruct
the scanned object from information about the locations of the X-ray
sources and measurements of their damping.

After discretizations (by the quadrature method or the Galerkin
method), we obtain a linear system of equations :math:`Ax=b`. The
generated test problems in Matrix Depot has type ``RegProb``, which is
defined as::

  immutable RegProb{T}
    A::Matrix{T}  # matrix of interest
    b::Vector{T}  # right-hand side
    x::Vector{T}  # the solution to Ax = b
  end

Here is an example::

  julia> matrixdepot("deriv2") # check information
  Computation of the Second Derivative:
             
  A classical test problem for regularization algorithms.
             
  Input options:
             
  [type,] n: the dimension of the matrix is n.
             
  Reference: P.C. Hansen, Regularization tools: A MATLAB pacakge for 
             analysis and solution of discrete ill-posed problems.

  
  julia> r = matrixdepot("deriv2", 3) # generate the test problem
  Test problems for Regularization Method
  A:
  3x3 Array{Float64,2}:
  -0.0277778   -0.0277778  -0.00925926
  -0.0277778   -0.0648148  -0.0277778 
  -0.00925926  -0.0277778  -0.0277778 
  b:
  [-0.01514653483985129,-0.03474793286789414,-0.022274315940957783]
  x:
  [0.09622504486493762,0.28867513459481287,0.48112522432468807]

  julia> r.A  # generate the matrix A
  3x3 Array{Float64,2}:
  -0.0277778   -0.0277778  -0.00925926
  -0.0277778   -0.0648148  -0.0277778 
  -0.00925926  -0.0277778  -0.0277778 

Here is a list test problems in the collection:

* :term:`deriv2`

.. glossary::
   :sorted:
      
   deriv2 
      Computation of the second derivative. The kernel :math:`K`
      is Green's function for the second derivative 

      .. math:: 

           K(s,t) = \begin{cases}
                    s(t - 1), \quad s < t, \\
                    t(s - 1), \quad s \geq t, \\
                    \end{cases}

      and both integration intervals are :math:`[0,1]`. The function 
      :math:`g` and :math:`f` are given by 

      .. math::

           g(s) = (s^3 - s)/6, \quad f(t) = t.

      The symmetric matrix :math:`A` and vectors :math:`x` and :math:`b` 
      are computed from :math:`K,f` and :math:`g` using the Galerkin method.
