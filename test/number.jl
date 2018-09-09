n = rand(1:38)
@test mdinfo(builtin(n)) != nothing
@test mdinfo(builtin(1:n)) != nothing
