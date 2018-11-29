.. _graphs:

Random Graphs
==============

* :term:`erdrey`
* :term:`gilbert`
* :term:`smallworld`

.. glossary::
   :sorted:

   erdrey
       An adjacency matrix of an Erdős–Rényi random graph: 
       an undirected graph is chosen uniformly at random from the set
       of all symmetric graphs with a fixed number of nodes and edges.
       For example::

          julia> using Random; Random.seed!(0);

          julia> matrixdepot("erdrey", Int8, 5, 3)
          5×5 SparseMatrixCSC{Int8,Int64} with 6 stored entries:
            [2, 1]  =  1
            [4, 1]  =  1
            [1, 2]  =  1
            [1, 4]  =  1
            [5, 4]  =  1
            [4, 5]  =  1

   gilbert
         An adjacency matrix of a Gilbert random graph: each possible edge occurs
	 independently with a given probability.

   smallworld
         Motivated by the small world model proposed by Watts and Strogatz [wast98]_, 
	 we proposed a random graph model by adding shortcuts to a kth nearest 
	 neighbor ring (node :math:`i` and :math:`j` are connected iff :math:`|i-j| \leq k` or 
	 :math:`|n - |i-j|| \leq k`).

.. code::

        julia> mdinfo("smallworld")
          Small World Network (smallworld)
          ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

          Generate an adjacency matrix for a small world network. We model it by adding shortcuts to a
          kth nearest neighbour ring network (nodes i and j are connected iff |i -j| <= k or |n - |i
          -j|| <= k.) with n nodes.

          Input options:

            •    [type,] n, k, p: the dimension of the matrix is n. The number of nearest-neighbours
                to connect is k. The probability of adding a shortcut in a given row is p.

            •    [type,] n: k = 2 and p = 0.1.

          References:

         .. [wast98] D.J. Watts and S. H. Strogatz. Collective Dynamics of Small World 

		   Networks, Nature 393 (1998), pp. 440-442. 

