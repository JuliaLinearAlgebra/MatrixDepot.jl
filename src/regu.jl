# Test matrices for regularization methods from Hansen's
# Regularization toolbox

# compute a oscillating matrix from a given signular value spectrum
# using the idea proposed by Per Christian Hansen, "Test matrices for
# regularization methods. 
function oscillate{T}(Σ::Vector{T})
    n = length(Σ)
    dv = rand(T, 1, n)[:] + eps(T)
    ev = rand(T, 1, n-1)[:] + eps(T)
    B = Bidiagonal(dv, ev, true)
    U, S, V = svd(B)
    return U*diagm(Σ)*U'
end
        

#  mode = 1: geometrically distributed singular values.
#  mode = 2: arithmetrically distributed singular values.
#  κ = sqrt(1/eps(T)) is the condition number of the matrix.
function oscillate{T}(::Type{T}, n::Int, mode::Int)
    κ = sqrt(1/eps(T))
    if mode == 1
        factor = κ^(-1/(n-1))
        Σ = factor.^[0:n-1;]
    elseif mode == 2
        Σ = ones(T, n) - T[0:n-1;]/(n-1)*(1 - 1/κ)
    else
        error("invalid mode value.")
    end
    return oscillate(Σ)
end
oscillate{T}(::Type{T}, n::Int) = oscillate(T, n, 2)


immutable RegProb{T}
    A::AbstractMatrix{T}  # matrix of interest
    b::AbstractVector{T}  # right-hand side
    x::AbstractVector{T}  # the solution to Ax = b
end

function show(io::IO, p::RegProb)
    println(io, "Test problems for Regularization Methods")
    println(io, "A:\n", p.A)
    println(io, "b:\n", p.b)
    println(io, "x:\n", p.x)
end

# The following test problems are derived from Per Christian Hansen's
# Regularization tools for MATLAB. 
# http://www.imm.dtu.dk/~pcha/Regutools/
#
# BSD License
# Copyright (c) 2015, Per Christian Hansen
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in
#       the documentation and/or other materials provided with the distribution
#     * Neither the name of the DTU Compute nor the names
#       of its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


#
# Computation of the Second Derivative 
# 
function deriv2{T}(::Type{T}, n::Int, example::Int, matrixonly::Bool = true)
    h = convert(T, one(T)/n); sqh = sqrt(h) 
    h32 = h*sqh; h2 = h^2; sqhi = one(T)/sqh
    t = 2/3; A = zeros(T, n, n)

    # compute A
    for i = 1:n
        A[i,i] = h2*((i^2 - i + 0.25)*h - (i - t))
        for j = 1:i-1
            A[i,j] = h2*(j - 0.5)*((i - 0.5)*h - 1)
        end
    end
    A = A + tril(A, -1)'
    
    if matrixonly
        return A
    else
        b = zeros(T, n)
        x = zeros(T, n)
        if example == 1
            # compute b
            [b[i] = h32*(i - 0.5)*((i^2 + (i-1)^2)*h2/2 - 1)/6 for i = 1:n]
        
            # compute x
            [x[i] = h32*(i - 0.5) for i = 1:n]
        elseif example == 2
            ee = one(T) - exp(one(T))
            [b[i] = sqhi*(exp(i*h) - exp((i-1)*h) + ee*(i-0.5)*h2 - h) for i = 1:n]
            [x[i] = sqhi*(exp(i*h) - exp((i-1)*h)) for i = 1:n]
        elseif example == 3
            mod(n, 2) == 0 || error("The order n must be even.")
            for i = 1:div(n,2)
                s12 = (i*h)^2; s22 = ((i-1)*h)^2
                b[i] = sqhi*(s12 + s22 - 1.5)*(s12 - s22)/24
            end
            for i = div(n, 2)+1:n
                s1 = i*h; s12 = s1^2; s2 = (i-1)*h; s22 = s2^2
                b[i] = sqhi*(-(s12 + s22)*(s12 - s22) + 4*(s1^3 - s2^3) 
                       - 4.5*(s12 - s22) + h)/24
            end
            [x[i] = sqhi*((i*h)^2 - ((i-1)*h)^2)/2 for i = 1:div(n, 2)]
            [x[i] = sqhi*(h - ((i*h)^2 - ((i-1)*h)^2)/2) for i = div(n, 2)+1:n]
        else
            error("Illegal value of example.")
        end       
        return RegProb(A, b, x)
    end
end

deriv2{T}(::Type{T}, n::Int, matrixonly::Bool = true) = deriv2(T, n, 1, matrixonly)


#
# One-Dimensional Image Restoration Model
# 
function shaw{T}(::Type{T}, n::Int, matrixonly::Bool = true)
    mod(n, 2) == 0 || error("The dimension of the matrix must be even.")
    h = pi/n; A = zeros(T, n, n)
    
    # compute A
    co = cos(-pi/2 + T[.5:n-.5;]*h)
    psi = pi*sin(-pi/2 + T[.5:n-.5;]*h)
    for i = 1:div(n,2)
        for j = i:n-i
            ss = psi[i] +psi[j]
            A[i,j] = ((co[i] + co[j])*sin(ss)/ss)^2
            A[n-j+1, n-i+1] = A[i,j]
        end
        A[i, n-i+1] = (2*co[i])^2
    end
    A = A + triu(A, 1)'; A = A*h
    
    if matrixonly
        return A
    else
        # compute x and b
        a1 = 2; c1 = 6; t1 = .8
        a2 = 1; c2 = 2; t2 = -.5
        x = a1*exp(-c1*(-pi/2 + T[.5:n-.5;]*h - t1).^2) + 
            a2*exp(-c2*(-pi/2 + T[.5:n-.5;]*h - t2).^2)
        b = A*x
        return RegProb(A, b, x)
    end
end

#
# A Problem with a Discontinuous Solution
#
function wing{T}(::Type{T}, n::Int, t1::Real, t2::Real, matrixonly = true)
    t1 < t2 || error("t1 must be smaller than t2")
    A = zeros(T, n, n); h = 1/n
    
    # compute A
    sti = (T[1:n;]-0.5)*h
    [A[i,:] = h*sti.*exp(-sti[i] * sti.^2) for i = 1:n]

    if matrixonly
        return A
    else
        # compute b
        b = sqrt(h)*0.5*(exp(-sti*t1^2) - exp(-sti*t2^2))./sti

        # compute x
        indices = [findfirst(t1 .< sti, true): findlast(t2 .> sti, true);]
        x = zeros(T,n); x[indices] = sqrt(h)*ones(length(indices))
        return RegProb(A, b, x)
    end
end
wing{T}(::Type{T}, n::Int, matrixonly = true) = wing(T, n, 1/3, 2/3, matrixonly)

#
# Severely Ill-posed Problem Suggested by Fox & Goodwin
#
function foxgood{T}(::Type{T}, n::Int, matrixonly = true)
    h = 1/n; t = h*(T[1:n;] - one(T)/2)

    A = h*sqrt((t.^2)*ones(T,n)' + ones(T, n) * (t.^2)')

    if matrixonly
        return A
    else
        x = t
        b = zeros(T, n)
        [b[i] = ((one(T) + t[i]^2)^1.5 - t[i]^3)/3 for i in 1:n]
        return RegProb(A, b, x)
    end
end

#
# Inverse Heat Equation
#
function heat{T}(::Type{T}, n::Int, κ::Real, matrixonly = true)
    mod(n, 2) == 0 || error("The dimension of the matrix must be even.")
    h = one(T)/n; t = T[h/2:h:1;]
    c = h/(2*κ*sqrt(pi))
    d = one(T)/(4*κ^2)

    # compute the matrix A
    m = length(t); k = zeros(T, m)
    [k[i] = c*t[i]^(-1.5)*exp(-d/t[i]) for i in 1:m]
    r = zeros(T, m); r[1] = k[1]
    A = toeplitz(T, k, r)

    if matrixonly
        return A
    else
        # compute the vectors x and b
        x = zeros(T, n)
        for i = 1:div(n,2)  
            ti = i*20/n
            if ti < 2
                x[i] = 0.75*ti^2/4
            elseif ti < 3
                x[i] = 0.75 + (ti - 2)*(3 - ti)
            else
                x[i] = 0.75*exp(-(ti - 3)*2)
            end
        end
        x[div(n,2)+1:n] = zeros(T, div(n, 2))
        b = A*x
        return RegProb(A, b, x)
    end
end
heat{T}(::Type{T}, n::Int, matrixonly = true) = heat(T, n, 1, matrixonly)

#
# First Kind Fredholm Integral Equation
#
function baart{T}(::Type{T}, n::Int, matrixonly = true)
    mod(n, 2) == 0 || error("The dimension of the matrix must be even.")
    hs = pi/(2*n); ht = pi/n; c = one(T)/(3*sqrt(2))
    ht = convert(T, ht)
    A = zeros(T, n, n); ihs = T[0:n;]*hs; n1 = n+1; nh = div(n,2)
    f3 = exp(ihs[2:n1]) - exp(ihs[1:n])
    
    # compute A
    for j = 1:n
        f1 = f3; co2 = cos((j - one(T)/2)*ht); co3 = cos(j*ht)
        f2 = (exp(ihs[2:n1]*co2) - exp(ihs[1:n]*co2))/co2
        j == nh ? f3 = hs*ones(T, n) : 
                  f3 = (exp(ihs[2:n1]*co3) - exp(ihs[1:n]*co3))/co3
        A[:,j] = c*(f1 + 4*f2 + f3)
    end
    
    if matrixonly
        return A
    else
        # compute vector b
        si = T[.5:.5:n;]*hs; si = sinh(si)./si
        b = zeros(T, n)
        b[1] = 1 + 4*si[1] + si[2]
        b[2:n] = si[2:2:2*n-2] + 4*si[3:2:2*n-1] + si[4:2:2*n]
        b = b*sqrt(hs)/3

        # compute vector x
        x = -diff(cos(T[0:n;]*ht))/sqrt(ht)
        return RegProb(A, b, x)
    end
end

#
# Phillips's "famous" problem
#
function phillips{T}(::Type{T}, n::Int, matrixonly = true)
    mod(n, 4) == 0 || error("The dimension of the matrix must be a multiple of 4.")

    # compute A
    h = 12/n; n4 = div(n, 4); r1 = zeros(T,n)
    c = cos(T[-1:n4;]*4*pi/n)
    [r1[i] = h + 9/(h*pi^2)*(2*c[i+1] - c[i] - c[i+2]) for i in 1:n4]
    r1[n4+1] = h/2 + 9/(h*pi^2)*(cos(4*pi/n)-1)
    A = toeplitz(T, r1)
    

    if matrixonly
        return A
    else
        # compute the vector b
        b = zeros(T, n); c = pi/3
        for i = div(n,2)+1:n
            t1 = -6 + i*h; t2 = t1 - h
            b[i] = t1*(6-abs(t1)/2) + 
                ((3-abs(t1)/2)*sin(c*t1) - 2/c*(cos(c*t1) - one(T)))/c -  
                t2*(6-abs(t2)/2) -
                ((3-abs(t2)/2)*sin(c*t2) - 2/c*(cos(c*t2) - one(T)))/c
            b[n-i+1] = b[i]
        end
        [b[i] = b[i]/sqrt(h) for i=1:n]
    
        # compute x
        x = zeros(T, n)
        x[2*n4+1:3*n4] = (h + diff(sin(T[0:h:(3+10*eps(T));]*c))/c)/sqrt(h)
        x[n4+1:2*n4] = x[3*n4:-1:2*n4+1]
        return RegProb(A, b, x)
    end
end

# replicates the grid vectors xgv and ygv to produce a full grid. 
function meshgrid(xgv, ygv)
    X = [i for j in ygv, i in xgv]
    Y = [j for j in ygv, i in xgv]
    return X, Y
end

#
# one-dimensional gravity surveying problem
#
function gravity{T}(::Type{T}, n::Int, example::Int, 
                    a::Number, b::Number, d::Number, matrixonly::Bool = true)
    dt = one(T)/n
    a = convert(T, a); b = convert(T, b); d = convert(T, d)
    ds = (b-a)/n
    tv = dt*(T[1:n;] - one(T)/2)
    sv = a + ds*(T[1:n;] - one(T)/2)
    Tm,Sm = meshgrid(tv,sv)
    A = dt*d*ones(T, n, n)./(d^2 + (Sm-Tm).^2).^(convert(T, 3/2))
    if matrixonly
        return A
    else
        x = ones(T, n)
        nt = @compat round(Integer, n/3)
        nn = @compat round(Integer, n*7/8)
        if example == 1
            x = sin(pi*tv) + 0.5*sin(2*pi*tv)
        elseif example == 2
            x[1:nt] = (2/nt)*[1:nt;]
            x[nt+1:nn] = ((2*nn - nt) - [nt+1:nn;])/ (nn-nt)
            x[nn+1:n] = (n - [nn+1:n;])/(n-nn)
        elseif example == 3
            x[1:nt] = 2*ones(T, nt)
        else
            error("Illegal value of example")
        end
        b = A*x
        return RegProb(A, b, x) 
    end
end

gravity{T}(::Type{T}, n::Int, example::Int, matrixonly::Bool = true) = 
    gravity(T, n, example, 0, 1, 0.25, matrixonly)

gravity{T}(::Type{T}, n::Int, matrixonly::Bool = true) = 
           gravity(T, n, 1, 0, 1, 0.25, matrixonly) 

#
# Image deblurring test problem
#
function blur{T}(::Type{T}, n::Int, band::Int, σ::Number, 
                 matrixonly::Bool = true)
    σ = convert(T, σ)
    z = [exp(-(T[0:band-1;].^2)/(2*σ^2)); zeros(T, n-band)]
    A = toeplitz(T, z)
    A = sparse(A)
    A = (1/(2*pi*σ^2))*kron(A, A)
    if matrixonly
        return A
    else
        # start with an image of all zeros
        if VERSION < v"0.4.0-dev+1419"
            n2 = @compat round(Integer, n/2)
            n3 = @compat round(Integer, n/3)
            n6 = @compat round(Integer, n/6)
            n12 = @compat round(Integer, n/12)
        else
            n2 = round(Integer, n/2, RoundNearestTiesAway)
            n3 = round(Integer, n/3, RoundNearestTiesAway)
            n6 = round(Integer, n/6, RoundNearestTiesAway)
            n12 = round(Integer, n/12, RoundNearestTiesAway)
        end
        m = max(n, 2*n6+1+n2+n12)
        x = zeros(T, m, m)
        
        # add a large ellipse
        Te = zeros(T, n6, n3)
        [if ((i/n6)^2 + (j/n3)^2 < 1) Te[i,j] = 1 end for i = 1:n6, j = 1:n3]
        Te = [flipdim(Te, 2) Te;]
        Te = [flipdim(Te, 1); Te]

        x[2+[1:2*n6;], n3-1+[1:2*n3;]] = Te

        
        # add a smaller ellipse
        Te = zeros(T, n6, n3)
        [if ((i/n6)^2 + (j/n3)^2 < 0.6) Te[i,j] = 1 end for i = 1:n6, j = 1:n3]
        Te = [flipdim(Te, 2) Te;]
        Te = [flipdim(Te, 1); Te]
        x[n6+[1:2*n6;], n3-1+[1:2*n3;]] = x[n6+[1:2*n6;],n3-1+[1:2*n3;]] + 2*Te
        f = find((i) -> i==3, x)
        x[f] = 2*ones(T, length(f))
        
        # Add a triangle
        Te = triu(ones(T, n3, n3))
        mT, nT = size(Te)
        x[n3+n12+[1:nT;], 1+[1:mT;]] = 3*Te


        # add a cross
        Te = zeros(T, 2*n6+1, 2*n6+1)
        mT, nT = size(Te)
        Te[n6+1, 1:nT] = ones(T, nT)
        Te[1:mT, n6+1] = ones(T, mT)
        x[n2+n12+[1:mT;], n2+[1:nT;]] = 4*Te

        x = reshape(x[1:n, 1:n], n^2)
        b = A*x 
        return RegProb(A, b, x) 
    end
end
blur{T}(::Type{T}, n::Int, matrixonly::Bool = true) = blur(T, n, 3, 0.7, matrixonly)


#
# Inverse Laplace transformation
#
function ilaplace{T}(::Type{T}, n::Int)

end

#
# Stellar parallax problem with real observations
#
function parallax{T}(::Type{T}, n::Int)

end

#
# Test problem with a "spiky" solution
#
function spikes{T}(::Type{T}, n::Int)
# need meshgrid
end

#
# Two-dimensional tomography problem with sparse matrix
#
function tomo{T}(::Type{T}, n::Int)

end

#
# Integral equation with no square integrable solution
#
function ursell{T}(::Type{T}, n::Int)

end
