# Groups

```@setup groups
using MatrixDepot
```

Groups categorize the problems in Matrix Depot. Group names are symbols, e.g.
`:symmetric`, and can be used as patterns anywhere a pattern is accepted.
[`listgroups`](@ref) returns all defined groups.

## Predefined groups

| group         | meaning                                                      |
|:--------------|:-------------------------------------------------------------|
| `:all`        | all problems in the collection                               |
| `:builtin`    | locally generated built-in problems                          |
| `:user`       | user-defined generators, see [User Generators](user.md)         |
| `:local`      | `:builtin` and `:user`                                       |
| `:eigen`      | part of the eigensystem is explicitly known                  |
| `:graph`      | adjacency matrix of a graph                                  |
| `:illcond`    | ill-conditioned for some parameter values                    |
| `:inverse`    | the inverse is known explicitly                              |
| `:orthogonal` | the matrix is orthogonal                                     |
| `:posdef`     | positive definite for some parameter values                  |
| `:random`     | the matrix has random entries                                |
| `:regprob`    | test problem for regularization methods                      |
| `:sparse`     | the matrix is sparse                                         |
| `:symmetric`  | symmetric for some parameter values                          |

```@example groups
listnames(:posdef)
```

## User-defined groups

[`setgroup!`](@ref) defines a new group from any pattern, not only from lists of names,
and [`deletegroup!`](@ref) removes it again:

```@example groups
setgroup!(:myfav, ["lehmer", "cauchy", "hilb"])
setgroup!(:spd, :symmetric & :posdef)
mdlist(:myfav)
```

```@example groups
deletegroup!(:myfav)
deletegroup!(:spd)
```
