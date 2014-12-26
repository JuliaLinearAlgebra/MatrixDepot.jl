
.. _matrices:

Matrices in the Collection 
--------------------------

 
* :term:`circul`
* :term:`frank`
* :term:`grcar`
* :term:`hilb`
* :term:`rosser`


.. glossary::
   :sorted:

   grcar
      Grcar matrix is a Toeplitz matrix with sensitive eigenvalues. The
      image below is a 200-by-200 Grcar matrix used in [nrt92]_.

      .. image:: images/grcar.png

      .. [nrt92] N.M. Nachtigal, L. Reichel and L.N. Trefethen, A hybrid
		 GMRES algorithm for nonsymmetric linear system, SIAM J. 
		 Matrix Anal. Appl., 13 (1992), pp. 796-825.

   frank
      Frank matrix is an upper Hessenberg matrix with determinant 1. 
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
