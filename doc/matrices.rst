
.. _matrices:

Matrices in the Collection 
--------------------------

* :term:`cauchy` 
* :term:`circul`
* :term:`dingdong`
* :term:`forsythe`
* :term:`frank`
* :term:`grcar`
* :term:`hadamard`
* :term:`hilb`
* :term:`invhilb`
* :term:`magic`
* :term:`moler`
* :term:`rosser`
* :term:`sampling`
* :term:`triw`

.. glossary::
   :sorted:

   pascal
      


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
