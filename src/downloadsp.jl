
# Download extra data specific form suite sparse at TAMU
#

const SS_SITE = "https://sparse.tamu.edu"
const SS_FILES = "/files/ssstats.csv"
const SS_INDEX = "/files/ss_index.mat"
const SS_GRAPHICS = "/files/" # to be followed by /$Group/$(Name)$ending
# where endings indicate graphics about the main matrix of different kinds
# ".png"|"_thumb.png"               - sparsity structure 
# "_graph.gif"|"_graph_thumb.gif"   - iconnectivity graph
# "_svd.png"                        - singular values plot

const SS_SVDDIR = "/svd/" # to be followed by "$Group/$(Name)_SVD{.mat|_piroband.mat}

# the files "$SS_SITE/$Group/$Name" contain html-formatted further information 
# the directories "/mat" "/MM" "/RB"
# contain data files in $Group/$Name.mat or $Group/Name.tar.gz" in different formats
#

# file "/statistics" contains html with legend of data contained in "ss_index"

# the SVD-based statistics data are only available on the "further information html file"
# The may be derived from the complete list of SVD values from "/SVD" if available
#

#=

Statistics computed for the SuiteSparse Matrix Collection

LastRevisionDate
    This is a single string kept in the UF_index MATLAB struct that states when the collection or the index was last modified.
DownloadTimeStamp
    The date and time the index that you last downloaded the index.
Group
    A cell array. Group{id} is the group name for the matrix whose serial number is 'id'. Each matrix has a unique id number in the range of 1 to the number of matrices in the collection. Once an id is assigned to a matrix, it never changes.
Name
    Name{id} is the name of the matrix (excluding the Group). Name{id} is not unique. The full name of a matrix should always be given as Group/Name.
nrows
    The number of rows in the matrix.
ncols
    The number of columns in the matrix.
nnz
    The number of numerically nonzero entries in the matrix, or nnz(Problem.A) in MATLAB, where Problem=UFget(id) is a struct containing the MATLAB format of the problem. This statistic can differ from the number of 'entries' explicitly stored in the matrix, however, since some of these entries may be numerically zero. In the MATLAB format, these explicit zero entries are stored in the binary Problem.Zeros matrix, since MATLAB drops all explicit zeros from its sparse matrix storage. The Problem.A matrix in MATLAB has nnz entries in it, with no explicit zeros. In the Matrix Market and Rutherford-Boeing format, a single file holds all entries, both nonzero and the explicit zero entries.
nzero
    The number of explicit entries present in the matrix that are provided by the matrix author but which are numerically zero. nzero(id) is nnz(Problem.Zeros).
pattern_symmetry
    Let S=spones(Problem.A) be the binary pattern of the explicit nonzeros in the matrix. Let pmatched be the number of matched offdiagonal entries, where both S(i,j) and S(j,i) are one, with i not equal to j. Let nzoffdiag be the number of offdiagonal entries in S. Then pattern_symmetry is the ratio of pmatched/nzoffdiag. Note that if S(i,j) and S(j,i) are both one, then this pair of entries is counted twice in both pmatched and nzoffdiag. If the matrix is rectangular, this statistic is zero. If there are no offdiagonal entries, the statistic is 1.
numerical_symmetry
    Let xmatched be the number of matched offdiagonal entries, where A(i,j) is equal to the complex conjugate of A(j,i) and where i and j are not equal. Then numerical_symmetry is the ration xmatched / nzoffdiag (or 1 if nzoffdiag is zero). This statistic is zero for rectangular matrices. Note that this statistic measures how close a matrix is to being Hermitian (A=A' in MATLAB). For complex symmetric matrices (A=A.' in MATLAB), this ratio will be less than one (unless all offdiagonal entries are real).
isBinary
    1 if the matrix is binary (all entries are 0 or 1), 0 otherwise.
isReal
    1 if the matrix is real, 0 if complex.
nnzdiag
    The number of numerically nonzero entries on the diagonal (nnz (diag (Problem.A)) in MATLAB notation). This statistic ignores explicit zero entries (Problem.Zeros in the MATLAB struct).
posdef
    1 if the matrix is known to be symmetric positive definite (or Hermitian positive definite for the complex case), 0 if it is known not to be, -1 if it is symmetric (or Hermitian) but with unknown positive-definiteness. If the statistic is unknown (-1) this may be revised in subsequent versions of the index.
amd_lnz
    Let C=S+S' where S=spones(A) is the binary pattern of A. Then amd_lnz is number of nonzeros in the Cholesky factorization of the matrix C(p,p) (assuming C is positive definite) where p=amd(C) is the AMD fill-reducing ordering. This statistic is -2 for rectangular matrices or for graph problems. It is -1 if it is not computed. This figure gives an estimate of the memory requirements for x=A\b in MATLAB for square matrices.
amd_flops
    The floating-point operation count for computing the Cholesky factorization of C(p,p) (see above).
amd_vnz
    The number of entries in a Householder-vector representation of the Q factor of R (but not the QR in MATLAB), after a COLAMD fill-reducing ordering. This is an upper bound on L for the LU factorization of A.
amd_rnz
    The number of entries in R for the QR factorization of A, after a COLAMD fill-reducing ordering. This is an upper bound on U for the LU factorization of A.
nblocks
    The number of blocks from dmperm (see dmperm in MATLAB).
sprank
    The structural rank of the matrix, which is the number of maximual number of nonzero entries that can be permuted to the diagonal (see dmperm, or sprank in MATLAB). This statistic is not computed for problems that represent graphs, since in those cases the diagonal of the matrix is often not relevant (self-edges are often ignored).
RBtype
    The Rutherford Boeing type of the matrix (ignoring explicit zeros in Problem.Zeros). RBtype is a a 3-letter lower-case string. The first letter is:

    r
        if A is real but not binary
    p
        if A is binary
    c
        if A is complex
    i
        if A is integer

    The second letter:

    r
        if A is rectangular
    u
        if A is square and unsymmetric
    s
        if A is symmetric (nnz(A-A.') is zero in MATLAB)
    h
        if A is Hermitian (nnz(A-A') is zero in MATLAB)
    z
        if A is skew-symmetric (nnz(A+A.') is zero in MATLAB)

    The third letter is always 'a' (for 'assembled'). The RB format allows for unassembled finite-element matrices, but they are converted to assembled format for this collection.
cholcand
    1 if the matrix is symmetric (Hermitian if complex) and if all the diagonal entries are postive and real. Zero otherwise. If 1, the matrix is a candidate for a Cholesky factorization, which MATLAB will first try when computing x=A\b.
ncc
    The number of of strongly-connected components in the graph of A (if A is square) or in the bipartite graph if A is rectangular. The diagonal is ignored.
isND
    1 if the problem comes from a 2D/3D discretization, zero otherwise. This determination is not a property of the matrix, but a qualitative assesment of the kind of problem the matrix represents.
isGraph
    1 if the problem is best considered as a graph rather than a system of equations, zero otherwise. This determination is not a property of the matrix, but a qualitative assesment of the kind of problem the matrix represents.

UFstats.csv

A CSV file is also available with some of this index information (UFstats.csv). The first line of the CSV file gives the number of matrices in the collection, and the second line gives the LastRevisionDate. Line k+2 in the file lists the following statistics for the matrix whose id number is k: Group, Name, nrows, ncols, nnz, isReal, isBinary, isND, posdef, pattern_symmetry, numerical_symmetry, and kind.
SVD-based statistics

The following statistics are not (yet) in the UFindex. They are currently available only on the web page for each matrix. You can also download the singular values at www.cise.ufl.edu/research/sparse/svd. These were typically calculated with s=svd(full(A)) in MATLAB, are are thus only available for modest-sized matrices.

norm(A)
    The 2-norm of A (the largest singular value)
min(svd(A))
    The smallest singular value
cond(A)
    The 2-norm condition number, which is the ratio of the largest over the smallest singular value.
rank(A)
    The rank of the matrix, which is the number of singular values larger than the tolerance of max(m,n)*eps(norm(A)). This tolerance is plotted in green in the figure.
sprank(A)-rank(A)
    sprank(A) (see above) is an upper bound on the rank of A.
null space dimension
    The dimension of the null space (zero if it has full numerical rank). This is simply min(m,n)-rank(A).
full numerical rank?
    'yes' or 'no'.
singular value gap
    If k=rank(A), and the matrix is rank deficient, then this is the ratio s(k)/s(k+1). A red line between the kth and (k+1)st singular value is drawn to illustrate this gap.
singular values
    These can be downloaded as a MATLAB MAT-file. Each file contains a struct with the fields: s (a vector containing the singular values), how (a string stating how the SVD was computed), and status (a string that is either 'ok' or a warning). If the status shows that the SVD did not converge, the singular values are probably not computed accurately.

=#

using MAT

const MFIELDNAMES = ["Group", "Name", "nrows", "ncols", "nnz", "nzero", "nnzdiag", 
                    "pattern_symmetry", "numerical_symmetry",
                    "isBinary", "isReal", "isND", "isGraph", "posdef", "cholcand",
                    "amd_lnz", "amd_vnz", "amd_rnz","amd_flops",
                    "nblocks", "ncc", "sprank", "RBtype",
                    "lowerbandwidth", "upperbandwidth",
                    "rcm_lowerbandwidth", "rcm_upperbandwidth",
                    "xmin", "xmax"]
const MFIELDTYPES = [String, String, Int, Int, Int, Int, Int,
                     Float64, Float64,
                     Bool, Bool, Bool, Bool, Bool, Bool,
                     Int, Int, Int, Int,
                     Int, Int, Int, String,
                     Int, Int, Int, Int,
                     ComplexF64, ComplexF64]

function read_ss_index()
    file = joinpath(DATA_DIR, "ss_index.mat")
    m = matopen(file) do io
        read(io, "ss_index")
    end
    d = Dict{String,Any}()
    for i in 1:length(MFIELDNAMES)
        fn = MFIELDNAMES[i]
        a = m[fn]
        T = MFIELDTYPES[i]
        b = vec([from_matdat(T, s) for s in a])
        push!(d, fn => b)
    end
    d
end

# Group, Name, nrows, ncols, nnz, isReal, isBinary, isND, posdef, pattern_symmetry, numerical_symmetry, and kind. 
const FIELDNAMES = ["Group", "Name", "nrows", "ncols", "nnzold", "isReal", "isBinary", "isND", "posdef", "pattern_symmetry", "numerical_symmetry", "kind", "nnz"] 
const FIELDTYPES = [String, String, Int, Int, Int, Bool, Bool, Bool, Bool, Float64, Float64, String, Int]

function read_ssstats()
    file = joinpath(DATA_DIR, "ssstats.csv")
    fieldno = length(FIELDNAMES)
    d = open(file) do io
        n = parse(Int, readline(io))
        d = makedict(n)
        a = [d[fn] for fn in FIELDNAMES]
        timestamp = readline(io)
        index = 0
        while (line = readline(io)) != ""
            fields = split(line, ',', limit=14)
            length(fields) >= fieldno || throw(DataError("invalid input line")) 
            index >= n && break
            index += 1
            for i in 1:fieldno
                s = fields[i]
                a[i][index] = from_string(FIELDTYPES[i], s)
            end
        end
        n == index || throw(DataError("missing data"))
        d
    end
    d
end

function makedict(n::Integer)
    d = Dict{String,Any}()
    for i in 1:length(FIELDNAMES)
        fn = FIELDNAMES[i]
        T = FIELDTYPES[i]
        push!(d, fn => Vector{T}(undef, n))
    end
    d
end

from_string(T::Type, s::AbstractString) = parse(T, s)
from_string(::Type{String}, s::AbstractString) = s
from_matdat(T::Type, s::Any) = T(s)
from_matdat(::Type{String}, s::AbstractString) = s

