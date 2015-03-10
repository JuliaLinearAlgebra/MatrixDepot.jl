# This is a modification of function mmread in MatrixMarket.jl 
# (https://github.com/JuliaSparse/MatrixMarket.jl)

function mmreader(mmfile::IO, infoonly::Bool=false)
#   Reads the contents of the Matrix Market file 'filename' into a matrix,
#   which will be either sparse or dense, depending on the Matrix Market format
#   indicated by 'coordinate' (coordinate sparse storage), or 'array' (dense
#   array storage).
#
#   If infoonly is true (default: false), only information on the size and
#   structure is returned from reading the header. The actual data for the
#   matrix elements are not parsed.
#    mmfile = open(filename,"r")

    #Read first line
    firstline = chomp(readline(mmfile))
    tokens = split(firstline)
    length(tokens)==5 || throw(ParseError(string("Not enough words on first line: ", firstline)))
    tokens[1]=="%%MatrixMarket" || throw(ParseError(string("Not a valid MatrixMarket header:", firstline)))
    (head1, rep, field, symm) = map(lowercase, tokens[2:5])
    head1=="matrix" || throw(ParseError("Unknown MatrixMarket data type: $head1 (only \"matrix\" is supported)"))
    eltype = field=="real"    ? Float64 :
             field=="complex" ? Complex128 :
             field=="pattern" ? Bool :
             throw(ParseError("Unsupported field $field (only real and complex are supported)"))
    

    #Skip all comments and empty lines
    ll   = readline(mmfile)
    while length(chomp(ll))==0 || (length(ll) > 0 && ll[1] == '%') ll = readline(mmfile) end

    #Read matrix dimensions (and number of entries) from first non-comment line
    dd     = int(split(ll))
    length(dd) >= (rep == "coordinate" ? 3 : 2) || throw(ParseError(string("Could not read in matrix dimensions from line: ", ll)))
    rows   = dd[1]
    cols   = dd[2]
    entries = rep == "coordinate" ? dd[3] : rows * cols
    if infoonly return rows, cols, entries, rep, field, symm end
    if rep == "coordinate"
        rr = Array(Int, entries)
        cc = Array(Int, entries)
        xx = Array(eltype, entries)
        for i in 1:entries
            flds = split(readline(mmfile))
            rr[i] = int(flds[1])
            cc[i] = int(flds[2])
            xx[i] = eltype==Complex128 ? Complex128(float64(flds[3]), float64(flds[4])) :
                    eltype==Float64 ? float64(flds[3]) :
                    true
        end
      
        return  sparse(rr, cc, xx, rows, cols)
        
    end
    reshape([float64(readline(mmfile)) for i in 1:entries], (rows,cols))
end


