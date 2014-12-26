
.. _matrices:

Matrices in the Collection 
--------------------------


* :term:`hilb`
* :term:`rosser`


.. glossary::
   :sorted:

 
   rosser 
      The Rosser matrix's eigenvalues are very close together so it is 
      a challenging matrix for many eigenvalue algorithms. 
      ``matrixdepot("rosser", 8, 2, 1)`` generates the test matrix used 
      in the paper [rlhk51]_. ``matrixdepot("rosser")`` are more general
      test matrices with similar property. 
      
      .. image:: images/rosser.png
   
      .. [rlhk51] Rosser, Lanczos, Hestenes and Karush, J. Res. Natl. 
		  Bur. Stand. Vol. 47 (1951), pp291-297. `Archive <https://archive.org/details/jresv47n4p291>`_	 

   hilb 
      The Hilbert matrix is a very ill conditioned matrix. But it is 
      symmetric positive definite and totally positive so it is not a good  
      test matrix for Gaussian elimination [high02]_ (Sec. 28.1).

      .. image:: images/hilb.png

      .. [high02] Nicholas J. Higham. Accuracy and Stability of 
		  Numerical Algorithms, SIAM, PA, USA. 2002.
      
