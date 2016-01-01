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
    A::AbstractMatrix{T}  # matrix of interest
    b::AbstractVector{T}  # right-hand side
    x::AbstractVector{T}  # the solution to Ax = b
  end

Here is an example::

  julia> matrixdepot("deriv2")
  Computation of the Second Derivative:
             
  A classical test problem for regularization algorithms.
  
  Input options:
             
  1. [type,] n, [matrixonly]: the dimension of the matrix is n. 
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.)
           
  Reference: P.C. Hansen, Regularization tools: A MATLAB pacakge for 
             analysis and solution of discrete ill-posed problems. 
             Numerical Algorithms, 6(1994), pp.1-35

  julia> A = matrixdepot("deriv2", 4) # generate the test matrix
  4x4 Array{Float64,2}:
  -0.0169271   -0.0195313  -0.0117188  -0.00390625
  -0.0195313   -0.0481771  -0.0351563  -0.0117188 
  -0.0117188   -0.0351563  -0.0481771  -0.0195313 
  -0.00390625  -0.0117188  -0.0195313  -0.0169271 

  julia> A = matrixdepot("deriv2", Float32, 4)
  4x4 Array{Float32,2}:
  -0.0169271   -0.0195313  -0.0117188  -0.00390625
  -0.0195313   -0.0481771  -0.0351563  -0.0117188 
  -0.0117188   -0.0351563  -0.0481771  -0.0195313 
  -0.00390625  -0.0117188  -0.0195313  -0.0169271 

  
  julia> r = matrixdepot("deriv2", 3, false) # generate the test problem
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
* :term:`blur`
* :term:`deriv2`
* :term:`foxgood`
* :term:`gravity`
* :term:`heat`
* :term:`parallax`
* :term:`phillips`
* :term:`shaw`
* :term:`spikes`
* :term:`ursell`
* :term:`wing`


.. glossary::
   :sorted:

   ursell
     Discretization of a Fredholm integral equation of the first kind
     with kernel `K` and right-hand side `g` given by

     .. math:: 

         K(s,t) = \frac{1}{s+t+1}, \quad g(s) = 1,
  
     where both integration intervals are :math:`[0,1]` [ursell]_. 

     .. [ursell] F. Ursell, Introduction to the theory of linear 
  		 integral equations, Chapter 1 in L. M. Delves and 
		 J. Walsh, Numerical Solution of Integral Equations, 
		 Clarendon Press, Oxford, 1974.

   spikes
     Artificially generated discrete ill-posed problem. 


   gravity 
     One-dimensional gravity surveying model
     problem. Discretization of a 1-D model problem in gravity
     surveying, in which a mass distribution f(t) is located at depth
     d, while the vertical component of the gravity field g(s) is
     measured at the surface. The resulting problem is a first-kind
     Fredholm integral equation with kernel
     
     .. math::

        K(s,t) = d(d^2 + (s-t)^2)^{-3/2}.


   blur
     Image deblurring test problem. It arises in connection with the 
     degradation of digital images by atmospheric turbulence blur, 
     modelled by a Gaussian point-spread function

     .. math::
	
	h(x,y) = \frac{1}{2\pi\sigma^2}\exp(-\frac{x^2+y^2}{2\sigma^2}).

     The matrix `A` is a symmetric :math:`n^2\times n^2` doubly block
     Toeplitz matrix, stored in sparse format.

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


   parallax 
      Stellar parallax problem with 26 fixed, real observations.
      The underlying problem is a Fredholm integral equation of the
      first kind with kernel

      .. math::

          K(s,t) = \frac{1}{\sigma \sqrt{2\pi}}\exp\Big(-\frac{1}{2}\big(\frac{s-t}{\sigma}\big)^2\Big),

      with :math:`\sigma = 0.014234` and it is discretized by means of
      a Galerkin method with n orthonormal basis functions. The
      right-hand side b consists of a measured distribution function of
      stellar parallaxes, and its length is fixed at 26; i.e., the
      matrix :math:`A` is :math:`26\times n`. The exact solution, which
      represents the true distribution of stellar parallaxes, is
      unknown.
