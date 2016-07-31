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

         julia> matrixdepot("erdrey", 5, 3) # an undirected graph with 5 nodes and 3 edges.
	 5x5 sparse matrix with 6 Float64 entries:
	   [3, 1]  =  1.0
	   [3, 2]  =  1.0
	   [1, 3]  =  1.0
	   [2, 3]  =  1.0
	   [5, 4]  =  1.0
	   [4, 5]  =  1.0

   gilbert
         An adjacency matrix of a Gilbert random graph: each possible edge occurs
	 independently with a given probability.

   smallworld
         Motivated by the small world model proposed by Watts and Strogatz [wast98]_, 
	 we proposed a random graph model by adding shortcuts to a kth nearest 
	 neighbor ring (node :math:`i` and :math:`j` are connected iff :math:`|i-j| \leq k` or 
	 :math:`|n - |i-j|| \leq k`).

         .. [wast98] D.J. Watts and S. H. Strogatz. Collective Dynamics of Small World 
		   Networks, Nature 393 (1998), pp. 440-442. 
