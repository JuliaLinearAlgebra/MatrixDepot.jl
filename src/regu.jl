# Test matrices for regularization methods from Hansen's
# Regularization toolbox
       

"""
Oscillating Matrix
==================
A matrix A is called oscillating if A is totally 
    nonnegative and if there exists an integer q > 0 such that 
    A^q is totally positive.

*Input options:*

+ [type,] Σ: the singular vaule spectrum of the matrix;

+ [type,] n, mode: n is the dimension of the matrix. 
        mode = 1: geometrically distributed singular values.
        mode = 2: arithmetrically distributed singular values.

+ [type,] n: mode = 1.

['symmetric','pos-def', 'random', 'eigen'] 

*Reference:* Per Christian Hansen, Test matrices for 
    regularization methods. SIAM J. SCI. COMPUT Vol 16, 
    No2, pp 506-512 (1995).
"""
function oscillate{T}(Σ::Vector{T})
    n = length(Σ)
    dv = rand(T, 1, n)[:] + eps(T)
    ev = rand(T, 1, n-1)[:] + eps(T)
    B = Bidiagonal(dv, ev, true)
    U, S, V = svd(B)
    return U*diagm(Σ)*U'
end
function oscillate{T}(::Type{T}, n::Int, mode::Int)
    κ = sqrt(1/eps(T))
    if mode == 1
        factor = κ^(-1/(n-1))
        Σ = factor.^[0:n-1;]
    elseif mode == 2
        Σ = ones(T, n) - T[0:n-1;]/(n-1)*(1 - 1/κ)
    else
        throw(ArgumentError("invalid mode value."))
    end
    return oscillate(Σ)
end
oscillate{T}(::Type{T}, n::Int) = oscillate(T, n, 2)
oscillate(args...) = oscillate(Float64, args...)

immutable RegProb{T}
    A::AbstractMatrix{T}  # matrix of interest
    b::AbstractVector{T}  # right-hand side
    x::AbstractVector{T}  # the solution to Ax = b
end

immutable RegProbNoSolution{T}
    A::AbstractMatrix{T}  # matrix of interest
    b::AbstractVector{T}  # right-hand side
end

function show(io::IO, p::RegProb)
    println(io, "Test problems for Regularization Methods")
    println(io, "A:")
    display(p.A)
    println(io, "b:")
    display(p.b)
    println(io, "x:")
    display(p.x)
end

function show(io::IO, p::RegProbNoSolution)
    println(io, "Test problems for Regularization Methods with No Solution")
    println(io, "A:")
    display(p.A)
    println(io, "b:")
    display(p.b)
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


"""
Computation of the Second Derivative
====================================
A classical test problem for regularization algorithms:
    This is a mildly ill-posed problem. It is a discretization
    of a first kind Fredholm integral equation whose kernel K
    is the Green's function for the second derivative.
 
*Input options:*

+ [type,] n, example, [matrixonly]: the dimension of the 
        matrix is n.  One choose between between the following right-hand
        g and solution f: 

 
        example = 1 gives g(s) = (s^3 - s)/6, f(t) = t.

        example = 2 gives g(s) = exp(s) + (1 -e)s - 1, f(t) = exp(t)

        example = 3 gives  g(s) = | (4s^3 - 3s)/24,  s < 0.5
                                  | (-4s^3 + 12s^2 - 9s + 1)/24, s>= 0.5
                           f(t) = | t, t < 0.5
                           f(t) = | 1- t, t >= 0.5.

        If `matrixonly = false`, the linear system A, b, x will be 
        generated.  (`matrixonly = true` by default.)             

+ [type,] n, [matrixonly]: example = 1.

*Reference:* P. C. Hansen, Regularization tools: A MATLAB pacakge for 
    analysis and solution of discrete ill-posed problems. 
    Numerical Algorithms, 6(1994), pp.1-35
"""
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
            throw(ArgumentError("Illegal value of example."))
        end       
        return RegProb(A, b, x)
    end
end
deriv2{T}(::Type{T}, n::Int, matrixonly::Bool = true) = deriv2(T, n, 1, matrixonly)
deriv2(args...) = deriv2(Float64, args...)


"""
One-Dimensional Image Restoration Model
=======================================
This test problem uses a first-kind Fredholm integral equation
    to model a one-dimentional image restoration situation.

*Input options:*

+ [type,] n, [matrixonly]: the dimesion of the matrix n must be even.
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.)

*Reference:* C. B. Shaw, Jr., Improvements of the resolution of 
    an instrument by numerical solution of an integral equation. 
    J. Math. Ana. Appl. 37 (1972), 83-112.
"""
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
shaw(args...) = shaw(Float64, args...)

"""
A Problem with a Discontinuous Solution
=======================================

*Input options:*

+ [type,] n, t1, t2, [matrixonly]: the dimension of matrix is n. 
    t1 and t2 are two real scalars such that 0 < t1 < t2 < 1. 
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.);

+ [type,] n, [matrixonly]: t1 = 1/3 and t2 = 2/3.

*Reference:* G. M. Wing, A Primer on Integral Equations of the 
    First Kind, SIAM, 1991, p. 109.
"""
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
wing(args...) = wing(Float64, args...)

"""
Severely Ill-posed Problem Suggested by Fox & Goodwin
=====================================================
This is a model problem discretized by simple quadrature, which does 
not satifiy the discrete Picard condition for the small singular values.

*Input options:*

+ [type,] n, [matrixonly]: the dimension of the matrix is n.
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.)

*Reference:* C. T. H. Baker, The Numerical Treatment of Integral
    Equations, Clarendon Press, Oxford, 1977, p. 665.
"""
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
foxgood(args...) = foxgood(Float64, args...)

"""
Inverse Heat Equation
=====================

*Input options:*

+ [type,] n, κ, [matrixonly]: the dimension of the matrix is n and n 
    must be even. κ controls the ill-conditioning of the matrix.
    (κ = 5 gives a well-conditioned problem and κ = 1 
     gives an ill conditoned problem). 
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.);

+ [type,] n, [matrixonly]: κ = 1.

*Reference:* A. S. Carasso, Determining surface temperatures 
    from interior observations, SIAM J. Appl. Math. 42 (1982), 558-574.
"""
function heat{T}(::Type{T}, n::Int, κ::Real, matrixonly::Bool = true)
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
heat{T}(::Type{T}, n::Int, matrixonly::Bool = true) = heat(T, n, 1, matrixonly)
heat(args...) = heat(Float64, args...)

"""
Fredholm Integral Equation of the Fisrt Kind
============================================

*Input options:*

+ [type,] n, [matrixonly]: the dimenstion of the matrix is n.
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.)

*Reference:* M. L. Baart, The use of auto-correlation for 
    pesudo-rank determination in noisy ill-conditioned linear-squares
    problems, IMA, J. Numer. Anal. 2 (1982), 241-247.
"""
function baart{T}(::Type{T}, n::Int, matrixonly::Bool = true)
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
baart(args...) = baart(Float64, args...)

"""
Phillips's \"famous\" problem
=============================

*Input options:*

+ [type,] n, [matrixonly]: the dimenstion of the matrix is n.
    If matrixonly = false, the linear system A, b, x will be generated. 
    (matrixonly = true by default.)

*Reference:* D. L. Phillips, A technique for the numerical 
    solution of certain integral equations of the first kind, J. ACM
    9 (1962), 84-97.
"""
function phillips{T}(::Type{T}, n::Int, matrixonly::Bool = true)
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
phillips(args...) = phillips(Float64, args...)

# replicates the grid vectors xgv and ygv to produce a full grid. 
function meshgrid(xgv, ygv)
    X = [i for j in ygv, i in xgv]
    Y = [j for j in ygv, i in xgv]
    return X, Y
end

# MATLAB rounding behavior 
# This is equivalent to RoundNearestTiesAway 
# and it can be used for both Julia v0.3, v0.4
function round_matlab(x::AbstractFloat)
    y = trunc(x)
    ifelse(x==y,y,trunc(2*x-y))
end
round_matlab{T<:Integer}(::Type{T}, x::AbstractFloat) = 
        @compat trunc(T,round_matlab(x))   


"""
One-dimensional gravity surverying problem
==========================================
Discretization of a 1-D model problem in gravity surveying, in 
    which a mass distribution f(t) is located at depth d, while the
    vertical component of the gravity field g(s) is measured at the
    surface. 

*Input options:*

+ [type,] n, example, a, b, d, [matrixonly]: n is the dimension
        of the matrix. Three examples are implemented.
        

       (a) example = 1 gives f(t) = sin(pi*t) + 0.5*sin(2*pi*t).

       (b) example = 2 gives f(t) = piecewise linear function.

       (c) example = 3 gives f(t) = piecewise constant function.
 
       The t integration interval is fixed to [0, 1], while the s 
       integration interval [a, b] can be specified by the user.
       The parameter d is the depth at which the magnetic deposit is 
       located. The larger the d, the faster the decay of the singular values. 

       If matrixonly = false, the linear system A, b, x will be generated.
       (matrixonly = true by default.)

+ [type,] n, example, [matrixonly]: a = 0, b = 1, d = 0.25;

+ [type,] n, [matrixonly]: example = 1, a = 0, b = 1, d = 0.25.            

*Reference:* G. M. Wing and J. D. Zahrt, A Primer on Integral 
            Equations of the First Kind, SIAM, Philadelphia, 1991, p. 17.
"""
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
        nt = round_matlab(Integer, n/3)
        nn = round_matlab(Integer, n*7/8)
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
gravity(args...) = gravity(Float64, args...)

"""
Image deblurring test problem
=============================
The generated matrix A is an n*n-by-n*n sparse, symmetric, 
           doubly block Toeplitz matrix that models blurring of an n-by-n 
           image by a Gaussian point spread function.

*Input options:*

+ [type,] n, band, σ, [matrixonly]: the dimension of the matrix
          is n^2. band is the half-bandwidth, only matrix elements within
          a distance band-1 from the diagonal are nonzero. σ controls the
          width of the Gaussin point spread function. The larger the σ, the 
          wider the function and the more ill posed the problem. 
          If matrixonly = false, the linear system A, b, x will be generated.
              (matrixonly = true by default.)

+ [type,] n, [matrixonly]: band = 3, σ = 0.7.
"""
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
        n2 = round_matlab(Integer, n/2)
        n3 = round_matlab(Integer, n/3)
        n6 = round_matlab(Integer, n/6)
        n12 = round_matlab(Integer, n/12)

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
blur(args...) = blur(Float64, args...)


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

"""
Test problem with \"spike\" solution
====================================
Artifically generated discrete ill-posed problem.

*Input options:*

+ [type,] n, t_max, [matrixonly]: the dimension of the 
              matrix is n. t_max controls the length of the pulse train.
              If matrixonly = false, the linear system A, b, x will be 
              generated. (matrixonly = true by default.) The solution x
              consists a unit step at t = .5 and a pulse train of spike
              of decreasing magnitude at t = .5, 1.5, 2.5, ...;

+ [type,] n, [matrixonly]: t_max = 5.
"""
function spikes{T}(::Type{T}, n::Int, t_max::Int, matrixonly::Bool = true)
    del = convert(T, t_max/n)
    
    # compute A
    t, sigma = meshgrid(T[del:del:t_max;], T[del:del:t_max;])
    A = sigma./(2*sqrt(pi*t.^3)).*exp(-(sigma.^2)./(4*t))
    
    if matrixonly
        return A
    else
        # compute b and x
        heights = 2*ones(T, t_max); heights[1] = 25
        heights[2] = 9; heights[3] = 5; heights[4] = 4; heights[5] = 3
        x = zeros(T, n); n_h = one(Integer)
        peak = convert(T, 0.5/t_max); peak_dist = one(T)/ t_max
        if peak < 1            
            n_peak = round_matlab(Integer, peak*n)
            x[n_peak] = heights[n_h]
            x[n_peak+1:n] = ones(T, n-n_peak)
            peak = peak + peak_dist; n_h = n_h + one(Integer)
        end
        while peak < 1
            n_peak = round_matlab(Integer, peak*n)
            x[n_peak] = heights[n_h]
            peak = peak + peak_dist; n_h = n_h + 1
        end
        b = A*x
        return RegProb(A, b, x) 
    end
end
spikes{T}(::Type{T}, n::Int, matrixonly::Bool = true) = 
         spikes(T, n, 5, matrixonly)
spikes(args...) = spikes(Float64, args...)

#
# Two-dimensional tomography problem with sparse matrix
#
function tomo{T}(::Type{T}, n::Int)

end

"""
Integral equation with no square integrable solution 
====================================================
Discretization of a first kind Fredholm integral equation with 
kernel `K` and right-hand side `g` given by
                     `K(s,t) = 1/(s+t+1), g(s) = 1`,
where both integration intervals are `[0, 1]`. The matrix `A`
is a Hankel matrix.

*Input options:*

+ [type,] n, [matrixonly]: the dimension of the matrix
              is n. If matrixonly = false, the right-hand side b will also
              be generated. (matrixonly = true by default).

*Reference:* F. Ursell, Introduction to the theory of linear
              integral equations., Chapter 1 in L. M. Delves & J. Walsh,
              Numerical Solution of Integral Equations, Clarendon Press, 
              1974.
"""
function ursell{T}(::Type{T}, n::Int, matrixonly::Bool = true)
    r = zeros(T, n); c = copy(r)
    for k = 1:n
        d1 = one(T) + (one(T) + k)/n
        d2 = one(T) + k/n
        d3 = one(T) + (k - one(T))/n
        c[k] = n*(d1*log(d1) + d3*log(d3) - 2*d2*log(d2))
        e1 = one(T) + (n+k)/n
        e2 = one(T) + (n+k-one(T))/n
        e3 = one(T) + (n+k-2)/n
        r[k] = n*(e1*log(e1) + e3*log(e3) - 2*e2*log(e2))
    end
    A = hankel(T, c, r)
    if matrixonly
        return A
    else
        # compute the right-hand side b
        b = ones(T,n)/convert(T, sqrt(n))
        return RegProbNoSolution(A, b) 
    end
end
ursell(args...) = ursell(Float64, args...)
