.. _regu:

Test Problems for Regularization Methods
========================================

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
method), we obtain a linear system of equations :math:`Ax=b`. 
All the regularization test problems are derived from 
discretizations of a Fredholm integral equation of the first kind.
Each generated test problem has type ``RegProb``, which is
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

Here is a list of test problems in the collection:

* :term:`baart`
* :term:`deriv2`
* :term:`foxgood`
* :term:`heat`
* :term:`phillips`
* :term:`shaw`
* :term:`wing`

.. glossary::
   :sorted:

   heat 
     Inverse heat equation [carasso82]_. It is a Volterra integral equation of
     the first kind with integration interval :math:`[0,1]`. The
     kernel :math:`K` is given by 

     .. math::

        K(s,t) = k(s-t),

     where 

     .. math::

	k(t) = \frac{t^{-3/2}}{2\kappa \sqrt{\pi}}\exp\big(-\frac{1}{4\kappa^2t}\big).

     :math:`\kappa` controls the ill-conditioning of the matrix :math:`A`. 
     :math:`\kappa = 1` (default) gives an ill-conditioned matrix and 
     :math:`\kappa = 5` gives a well-conditioned matrix. 

     .. [carasso82] A.S. Carasso, Determining surface temperatures
		    from interior observations,
		    SIAM J. Appl. Math. 42 (1982), 558-574.

   baart 
     Discretization of an artificial Fredholm integral equation of
     the first kind [baart82]_. The kernel :math:`K` is given by 

     .. math::

        K(s,t) = \exp(s \cos (t)).

     The right-hand side :math:`g` and the solution :math:`f` are given by 

     .. math:: 

        g(s)=2\frac{\sin (s)}{s}, \quad f(t) = \sin(t). 

     .. [baart82] M.L. Baart, The use of auto-correlation for pesudo-rank
		  determination in noisy ill-conditioned linear least-squares
		  problems, IMA, J. Numer. Anal. 2 (1982), 241-247.

   phillips
     Phillips's "famous" problem. Discretization of the "famous" Fredholm
     integral equation of the first kind deviced by D.L. Phillips [phillips62]_. 
     The kernel :math:`K` and solution :math:`f` are given by

     .. math::

	K(s,t) = \theta(s-t), \quad f(t) = \theta(t),

     where 

     .. math::

	\theta(x) = \begin{cases}
                     1+\cos(\frac{\pi x}{3}), & |x| < 3, \\
                     0,            & |x| \geq 3. \\
                    \end{cases}

     The right-hand side :math:`g` is given by

     .. math::

        g(s) = (6 - |s|)\Big( 1 + \frac{1}{2}\cos\big(\frac{\pi s}{3}\big)\Big) + \frac{9}{2 \pi}\sin\Big(\frac{\pi |s|}{3}\Big).

     Both integration intervals are :math:`[-6,6]`. 

     .. [phillips62] D.L. Phillips, A technique for the numerical solution 
		     of certain integral equations of the first kind, J. ACM
		     9 (1962), 84-97.


   foxgood
      A severely ill-posed problem suggested by Fox & Goodwin. This 
      is a model problem which does not satisfy the discrete Picard
      condition for the small singular values [baker77]_.

      .. [baker77] C.T.H Baker, The Numerical Treatment of Integral 
		   Equations, Clarendon Press, Oxford, 1977, p. 665.

   wing
      A problem with a discontinuous solution. The kernel :math:`K` is 
      given by 

      .. math::

	 K(s,t) = t \exp(-st^2),

      with both integration intervals are :math:`[0,1]`.  The functions 
      :math:`f` and :math:`g` are given as 

      .. math::

	 f(t) = \begin{cases} 
	         1, \quad t_1 < t < t_2, \\
		 0, \quad \mbox{otherwise},\\
                \end{cases}
	 \quad
	 g(s) = \frac{\exp(-st_1^2) - \exp(-st_2^2)}{2s}.

      Here :math:`0 < t_1 < t_2 < 1`. The matrix :math:`A` and two
      vectors :math:`x` and :math:`b` are obtained by Galerkin discretization
      with orthonormal basis functions defined on a uniform mesh.

   shaw
      One-dimensional image restoration model. This test problem 
      uses a first-kind Fredholm integral equation to model a one-dimensional
      image restoration situation. The kernel :math:`K` is given by

      .. math::

           K(s,t) = (\cos(s)+\cos(t))^2\big(\frac{\sin(u)}{u}\big)^2,

      where 

      .. math:: 

           u = \pi(\sin(s) + \sin(t)).

      Both integration intervals are :math:`[-\pi/2, \pi/2]`. The solution
      :math:`f` is given by

      .. math::

           f(t) = a_1 \exp(-c_1(t-t_1)^2) + a_2 \exp(-c_2(t-t_2)^2).

      :math:`K` and :math:`f` are discretized by simple quadrature to 
      produce the matrix :math:`A` and the solution vector :math:`x`. 
      The right-hand :math:`b` is computed by :math:`b=Ax`.

      
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
