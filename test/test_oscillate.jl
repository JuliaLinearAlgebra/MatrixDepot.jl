let mode, n, A, p, k, A
# test two properties of oscillating matrix
# For an n x n oscillating matrix A, the following holds:
# 1. A has n distinct and postive eigenvalues λ_1 > λ_2 > ... > λ_n > 0
# 2. The ith eigenvector, corresponding to λ_i in the above ordering has 
#    exactly i-1 sign changes. 

n = 8 # rand(2:10)
@test_throws ArgumentError matrixdepot("oscillate", n, 4)
for mode = 1:2
A = matrixdepot("oscillate", n, mode)

eva, evc = eigen(A)
@test all([ei > 0 for ei in eva])

p = sortperm(eva, rev = true)
evc = evc[:,p]

# compute the number of sign change
function num_sign_change(v)
    signv = sign.(v)
    change = signv[1]
    num = 0
    for i in signv
        if i != change
            num += 1
            change = i
        end
    end
    return num
end

k = rand(1:n)
@test num_sign_change(evc[:,k]) == k - 1
end

@test matrixdepot("oscillate", n) !== nothing
end
println("'oscillate' passed test...")
