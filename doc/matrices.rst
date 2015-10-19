
.. _matrices:

Matrices  
========

* :term:`binomial`
* :term:`cauchy`
* :term:`chebspec` 
* :term:`chow`
* :term:`circul`
* :term:`clement`
* :term:`dingdong`
* :term:`fiedler`
* :term:`forsythe`
* :term:`frank`
* :term:`grcar`
* :term:`hadamard`
* :term:`hankel`
* :term:`hilb`
* :term:`invhilb`
* :term:`invol`
* :term:`kahan`
* :term:`kms`
* :term:`lehmer`
* :term:`lotkin`
* :term:`magic`
* :term:`minij`
* :term:`moler`
* :term:`neumann`
* :term:`oscillate`
* :term:`parter`
* :term:`pascal`
* :term:`pei`
* :term:`poisson`
* :term:`prolate`
* :term:`randcorr`
* :term:`rando`
* :term:`randsvd`
* :term:`rohess`
* :term:`rosser`
* :term:`sampling`
* :term:`toeplitz`
* :term:`tridiag`
* :term:`triw`
* :term:`vand`
* :term:`wathen`
* :term:`wilkinson`


.. glossary::
   :sorted:

   hankel
     `Hankel matrix <https://en.wikipedia.org/wiki/Hankel_matrix>`_ is a 
     a matrix that is symmetric and constant across the anti-diagonals.
     For example::

       julia> matrixdepot("hankel", [1,2,3,4], [7,8,9,10])
       4x4 Array{Float64,2}:
       1.0  2.0  3.0   4.0
       2.0  3.0  4.0   8.0
       3.0  4.0  8.0   9.0
       4.0  8.0  9.0  10.0

   toeplitz
     `Toeplitz matrix <https://en.wikipedia.org/wiki/Toeplitz_matrix>`_ is 
     a matrix in which each descending diagonal from left to right 
     is constant. For example::

       julia> matrixdepot("toeplitz", [1,2,3,4], [1,4,5,6])
       4x4 Array{Int64,2}:
       1  4  5  6
       2  1  4  5
       3  2  1  4
       4  3  2  1

       julia> matrixdepot("toeplitz", [1,2,3,4])
       4x4 Array{Int64,2}:
       1  2  3  4
       2  1  2  3
       3  2  1  2
       4  3  2  1

 


   prolate
      A prolate matrix is a symmetric ill-conditioned Toeplitz matrix

      .. math::

	  A = \begin{bmatrix}
	      a_0 & a_1 & \cdots \\
              a_1 & a_0 & \cdots \\
              \vdots & \vdots & \ddots \\
              \end{bmatrix}

      such that :math:`a_0= 2w` and :math:`a_k = (\sin 2 \pi wk)/\pi k` for 
      :math:`k=1,2, \ldots` and :math:`0<w<1/2` [varah93]_.

      .. [varah93] J.M. Varah. The Prolate Matrix. Linear Algebra and Appl.
                  187:267--278, 1993.

   oscillate 
      A matrix :math:`A` is called oscillating if :math:`A` is
      totally nonnegative and if there exists an integer q > 0 such 
      that A^q is totally positive. An :math:`n \times n` oscillating 
      matrix :math:`A` satisfies:
 
      1. :math:`A` has :math:`n` distinct and positive eigenvalues
	 :math:`\lambda_1 > \lambda_2 > \cdots > \lambda_n > 0`. 
      2. The :math:`i` th eigenvector,  corresponding to :math:`\lambda_i`
         in the above ordering, has exactly :math:`i -1` sign changes. 

      This function generates a symmetric oscillating matrix, which is useful 
      for testing numerical regularization methods [hansen95]_. For example::
	
	julia> A = matrixdepot("oscillate", 3)
	3x3 Array{Float64,2}:
	0.98694    0.112794   0.0128399 
	0.112794   0.0130088  0.0014935 
	0.0128399  0.0014935  0.00017282

	julia> eig(A)
	([1.4901161192617526e-8,0.00012207031249997533,0.9999999999999983],
	3x3 Array{Float64,2}:
	0.0119607   0.113658  -0.993448 
	-0.215799   -0.969813  -0.113552 
	0.976365   -0.215743  -0.0129276)

      .. [hansen95] Per Christian Hansen, Test matrices for
                    regularization methods. SIAM J. SCI. COMPUT Vol 16, No2,
                    pp 506-512 (1995)

   wathen 
      Wathen Matrix is a sparse, symmetric positive, random matrix arose 
      from the finite element method [wath87]_. The generated matrix 
      is the consistent mass matrix for a regular
      `nx-by-ny` grid of 8-nodes.

      .. image:: images/wathen.png

      .. [wath87] A.J. Wathen, Realistic eigenvalue bounds for the Galerkin
		  mass matrix, IMA J. Numer. Anal., 7 (1987), pp. 449-457.

   kms 
      Kac-Murdock-Szego Toeplitz matrix [tren89]_.

      .. image:: images/kms.png

      .. [tren89] W.F. Trench, Numerical solution of the eigenvalue
                  problem for Hermitian Toeplitz matrices,
                  SIAM J. Matrix Analysis and Appl., 10 (1989),
                  pp. 135-146 

   rohess
      A random orthogonal upper Hessenberg matrix. The matrix is 
      constructed via a product of Givens rotations.

      .. image:: images/rohess.png

   randsvd
      Random matrix with pre-assigned singular values. See [high02]_ (Sec. 28.3).

      .. image:: images/randsvd.png

   rando
      A random matrix with entries -1, 0 or 1.

      .. image:: images/rando.png

   wilkinson
      The Wilkinson matrix is a symmetric tridiagonal matrix with pairs
      of nearly equal eigenvalues. The most frequently used case is 
      ``matrixdepot("wilkinson", 21)``.

      .. image:: images/wilkinson.png


   neumann
      A singular matrix from the discrete Neumann problem. This matrix
      is sparse and the null space is formed by a vector of ones [plem76]_.

      .. image:: images/neumann.png

      .. [plem76] R.J. Plemmons, Regular splittings and the discrete Neumann
		  problem, Numer. Math., 25 (1976), pp. 153-161.


   poisson
      A block tridiagonal matrix from Poisson's equation. This matrix is
      sparse, symmetric positive definite and has known eigenvalues. 

      .. image:: images/poisson.png

   randcorr
      A random correlation matrix is a symmetric positive semidefinite 
      matrix with 1s on the diagonal.

      .. image:: images/randcorr.png


   chow
      The Chow matrix is a singular Toeplitz lower Hessenberg matrix. The
      eigenvalues are known explicitly [chow69]_.

      .. image:: images/chow.png
 
      .. [chow69] T.S. Chow, A class of Hessenberg matrices with known 
		  eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.

   parter
      The Parter matrix is a Toeplitz and Cauchy matrix with singular values near
      :math:`\pi` [part86]_. 

      .. image:: images/parter.png

      .. [part86] S. V. Parter, On the distribution of the singular values
		  of Toeplitz matrices, Linear Algebra and Appl., 80 (1986),
		  pp. 115-130.

   lehmer
      The Lehmer matrix is a symmetric positive definite matrix. It is 
      totally nonnegative. The inverse is tridiagonal and explicitly known 
      [neto58]_.

      .. image:: images/lehmer.png
     
      .. [neto58] M. Newman and J. Todd, The evaluation of matrix inversion
		  programs, J. Soc. Indust. Appl. Math, 6 (1958), pp. 466-476.


   tridiag
      A group of tridiagonal matrices. ``matrixdepot("tridiagonal", n)``
      generate a tridiagonal matrix with 1 on the diagonal and -2 on the 
      upper- lower- diagonal, which is a symmetric positive definite 
      M-matrix. This matrix is also known as Strang's matrix, named after
      Gilbert Strang.

      .. image:: images/tridiag.png

   binomial
      A binomial matrix that arose from the example in [bmsz01]_.
      The matrix is a multiple of involutory matrix.

      .. image:: images/binomial.png

      .. [bmsz01] G. Boyd, C.A. Micchelli, G. Strang and D.X. Zhou,
		  Binomial matrices, Adv. in Comput. Math., 14 (2001), pp 379-391.
 

   minij
      A matrix with :math:`(i,j)` entry ``min(i,j)``. It is a symmetric
      positive definite matrix. The eigenvalues and eigenvectors are 
      known explicitly. Its inverse is tridiagonal. 

      .. image:: images/minij.png


   clement
      The Clement matrix [clem59]_ is a Tridiagonal matrix with zero diagonal entries.
      If ``k = 1``, the matrix is symmetric.
      
      .. image:: images/clement.png
		 
      .. [clem59] P.A. Clement, A class of triple-diagonal matrices for test
		  purposes, SIAM Review, 1 (1959), pp. 50-52.

   fiedler
      The Fiedler matrix is symmetric matrix with a dominant positive eigenvalue
      and all the other eigenvalues are negative. For explicit formulas for 
      the inverse and determinant, see [todd77]_.

      .. image:: images/fiedler.png

      .. [todd77] J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
		  Birkhauser, Basel, and Academic Press, New York, 1977, pp. 159.


   lotkin
      The Lotkin matrix is the Hilbert matrix with its first row altered
      to all ones. It is unsymmetric, ill-conditioned and has many negative
      eigenvalues of small magnitude [lotk55]_.

      .. image:: images/lotkin.png
      
      .. [lotk55] M. Lotkin, A set of test matrices, MTAC, 9, (1955), pp. 153-161.


   chebspec 
      Chebyshev spectral differentiation matrix. 
      If ``k = 0``,the generated matrix is nilpotent and a vector with 
      all one entries is a null vector. If ``k = 1``, the generated matrix
      is nonsingular and well-conditioned. Its eigenvalues have negative 
      real parts.

      .. image:: images/chebspec.png

   invol
      An involutory matrix, i.e., a matrix that is its own inverse. 
      See [hoca63]_.

      .. image:: images/invol.png

      .. [hoca63] A.S. Householder and J.A. Carpenter, The singular 
		  values of involutory and idempotent matrices, Numer. 
		  Math. 5 (1963), pp. 234-237.
      

   vand
      The Vandermonde matrix is defined in terms of scalars 
      :math:`\alpha_0, \alpha_1, \ldots, \alpha_n` by 

      .. math::

	      V(\alpha_0, \ldots, \alpha_n) = \begin{bmatrix}
                                               1 & 1 & \cdots & 1 \\
					       \alpha_0 & \alpha_1 & \cdots & \alpha_n \\
					       \vdots   & \vdots   &        & \vdots   \\
					       \alpha_0^n  & \alpha_1^n & \cdots & \alpha_n^n \\
			                       \end{bmatrix}.
       
      The inverse and determinant are known explicitly [high02]_. 

      .. image:: images/vand.png

   pei
      The Pei matrix is a symmetric matrix with known inversison [pei62]_. 

      .. image:: images/pei.png

      .. [pei62] M.L. Pei, A test matrix for inversion procedures, Comm. ACM, 5 (1962), pp. 508.
		 
   kahan
      The Kahan matrix is a upper trapezoidal matrix, i.e., the 
      :math:`(i,j)` element is equal to 0 if :math:`i > j`. The useful
      range of ``theta`` is :math:`0 < theta < \pi`. The diagonal is 
      perturbed by ``pert*eps()*diagm([n:-1:1])``.


      .. image:: images/kahan.png


   pascal 
      The Pascal matrix's anti-diagonals form the Pascal's
      triangle:: 
      
        julia> matrixdepot("pascal", 6)
	6x6 Array{Int64,2}:
	1  1   1   1    1    1
	1  2   3   4    5    6
	1  3   6  10   15   21
	1  4  10  20   35   56
	1  5  15  35   70  126
	1  6  21  56  126  252
 
      See [high02]_ (28.4).

	    
      .. image:: images/pascal.png


   sampling
      Matrices with application in sampling theory. A n-by-n nonsymmetric matrix
      with eigenvalues :math:`0, 1, 2, \ldots, n-1` [botr07]_. 

      .. image:: images/sampling.png

      .. [botr07]  L. Bondesson and I. Traat, A Nonsymmetric Matrix with Integer
		   Eigenvalues, Linear and Multilinear Algebra, 55(3)(2007), pp. 239-247.

   moler
      The Moler matrix is a symmetric positive definite matrix. It has
      one small eigenvalue.

      .. image:: images/moler.png


   triw
      Upper triangular matrices discussed by Wilkinson and others [gowi76]_.

      .. image:: images/triw.png

      .. [gowi76] G.H. Golub and J.H. Wilkinson, Ill-conditioned eigensystems
		  and the computation of the Jordan canonical form, SIAM Review,
		  18(4), (1976), pp. 578-619.
      

   forsythe
      The Forsythe matrix is a n-by-n perturbed Jordan block. 

      .. image:: images/forsythe.png


   cauchy
      The Cauchy matrix is an m-by-n matrix with :math:`(i,j)` element
      
      .. math::
	 
	 \frac{1}{x_i - y_i}, \quad x_i - y_i \ne 0,

      where :math:`x_i` and :math:`y_i` are elements of vectors :math:`x` 
      and :math:`y`.
      
      .. image:: images/cauchy.png

   magic
      The magic matrix is a matrix with integer entries such that the 
      row elements, column elements, diagonal elements and anti-diagonal 
      elements all add up to the same number. 

      .. image:: images/magic.png

   hadamard
      The Hadamard matrix is a square matrix whose entries are 1 or -1. It 
      was named after Jacques Hadamard. The rows of a Hadamard matrix 
      are orthogonal. 
      
      .. image:: images/hadamard.png

   dingdong
      The Dingdong matrix is symmetric Hankel matrix invented by Dr. F. N. Ris
      of IBM, Thomas J Watson Research Centre. The eigenvalues cluster 
      around :math:`\pi/2` and :math:`-\pi/2` [nash90]_. 

      .. image:: images/dingdong.png

      .. [nash90] J.C. Nash, Compact Numerical Methods for Computers: Linear
		  Algebra and Function Minimisation, second edition, Adam Hilger, 
		  Bristol, 1990 (Appendix 1).

   invhilb
      Inverse of the Hilbert Matrix.

      .. image:: images/invhilb.png

   grcar
      The Grcar matrix is a Toeplitz matrix with sensitive eigenvalues. The
      image below is a 200-by-200 Grcar matrix used in [nrt92]_.

      .. image:: images/grcar.png

      .. [nrt92] N.M. Nachtigal, L. Reichel and L.N. Trefethen, A hybrid
		 GMRES algorithm for nonsymmetric linear system, SIAM J. 
		 Matrix Anal. Appl., 13 (1992), pp. 796-825.

   frank
      The Frank matrix is an upper Hessenberg matrix with determinant 1. 
      The eigenvalues are real, positive and very ill conditioned [vara86]_.  

      .. image:: images/frank.png

      .. [vara86] J.M. Varah, A generalization of the Frank matrix, SIAM J. Sci. Stat. 
		  Comput., 7 (1986), pp. 835-839.
   

   circul
      A circulant matrix has the property that each row is obtained by
      cyclically permuting the entries of the previous row one step 
      forward.

      .. image:: images/circul.png
 
   rosser 
      The Rosser matrix's eigenvalues are very close together so it is 
      a challenging matrix for many eigenvalue algorithms. 
      ``matrixdepot("rosser", 8, 2, 1)`` generates the test matrix used 
      in the paper [rlhk51]_. ``matrixdepot("rosser")`` are more general
      test matrices with similar property. 
      
      .. image:: images/rosser.png
   
      .. [rlhk51] Rosser, Lanczos, Hestenes and Karush, J. Res. Natl. 
		  Bur. Stand. Vol. 47 (1951), pp. 291-297. `Archive <https://archive.org/details/jresv47n4p291>`_	 

   hilb 
      The Hilbert matrix is a very ill conditioned matrix. But it is 
      symmetric positive definite and totally positive so it is not a good  
      test matrix for Gaussian elimination [high02]_ (Sec. 28.1).

      .. image:: images/hilb.png

      .. [high02] Nicholas J. Higham. Accuracy and Stability of 
		  Numerical Algorithms, SIAM, PA, USA. 2002.
      



.. note:: 
   The images are generated using `Winston.jl <https://github.com/nolta/Winston.jl>`_ 
   's ``imagesc`` function.
