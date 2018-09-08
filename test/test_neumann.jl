n = 10 # rand(2:10)
@test matrixdepot("neumann", 1) == reshape([4.0], 1, 1)
@test matrixdepot("neumann", n) == matrixdepot("neumann", Float64, n)

n_mesh, n = n, n^2

B = zeros(Float64, n, n)

i = 0

for i_mesh in 1:n_mesh, j_mesh in 1:n_mesh
    global i = i + 1
    if i_mesh > 1
        j = i - n_mesh
    else
        j = i + n_mesh
    end

    B[i,j] = B[i,j] - 1
    
    if j_mesh > 1
        j = i - 1
    else
        j = i + 1
    end

    B[i,j] = B[i,j] - 1

    j = i
    B[i,j] = 4
    
    if j_mesh < n_mesh
        j = i + 1
    else
        j = i - 1
    end

    B[i,j] = B[i,j] - 1

    if i_mesh < n_mesh
        j = i + n_mesh
    else
        j = i - n_mesh
    end

    B[i,j] = B[i,j] - 1
end
    
A = matrixdepot("neumann", n_mesh)

@test Matrix(A) ≈ B

println("'neumann' passed test...")
