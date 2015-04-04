# return info comment strings for UF sparse matrix
function ufinfo(filename)
    mmfile = open(filename,"r")
    info = "\n"
    ll = readline(mmfile)
    while length(ll) > 0 && ll[1] == '%'       
        info = string(info, ll)
        ll = readline(mmfile)
    end
    return info
end
