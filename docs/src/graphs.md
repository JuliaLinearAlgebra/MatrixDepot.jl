# Random Graphs

The group `:graph` contains generators of adjacency matrices of random graphs. They are
returned as sparse matrices:

```@example
using MatrixDepot, Random
Random.seed!(0)
matrixdepot("erdrey", Int8, 5, 3)
```

## erdrey

```@docs
MatrixDepot.erdrey
```

## gilbert

```@docs
MatrixDepot.gilbert
```

## smallworld

```@docs
MatrixDepot.smallworld
```

