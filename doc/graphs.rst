.. _graphs:

Random Graphs
==============

* :term:`erdrey`

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
