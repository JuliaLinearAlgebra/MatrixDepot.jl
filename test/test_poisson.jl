n = rand(1:n)

@test matrixdepot("poisson", Float64, n) == matrixdepot("poisson", n)

A = Matrix(matrixdepot("poisson", n))

n_mesh, n = n, n^2


B = zeros(Float64, n, n)

i = 0
for i_mesh in 1 : n_mesh
    for j_mesh in 1 : n_mesh
        global i = i + 1
        
        if i_mesh > 1 
            j = i - n_mesh
            B[i,j] = -1
        end

        if j_mesh > 1
            j = i - 1
            B[i,j] = -1
        end
        
        j = i
        B[i,j] = 4

        if j_mesh < n_mesh
            j = i + 1
            B[i,j] = -1
        end

        if i_mesh < n_mesh
            j = i + n_mesh
            B[i,j] = -1
        end
    end
end

@test A ≈ B

println("'poisson' passed test...")

