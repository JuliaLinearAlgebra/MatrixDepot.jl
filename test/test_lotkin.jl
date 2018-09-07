n = 10 # rand(1:10)
A = matrixdepot("hilb", n)
B = matrixdepot("lotkin", n)

@test A[2:n, :] == B[2:n, :] 
@test B[1,:][:] == ones(n)
println("'lotkin' passed test...")
