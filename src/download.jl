# Download data from the University of Florida Sparse Matrix Collection

function download_html()
    isfile("matrices.html") || download("http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html", "matrices.html")
    matrixdata = {}
    open("matrices.html") do f
        for line in readlines(f)
            if contains(line, "MAT</a>")
               collection_name, matrix_name   # split the data
               push!(matrixdata, (collection_name, matrix_name)) 
            end
        end
    end
    return matrixdata
end

function download_matrix(name)
    collection_name, matrix_name = split('/')
    if (collection_name, matrix_name) in matrixdata
        fn = string(, ".mat") # combine to string
        if !isfile(fn)
            url = "$fn"
            
            try 
                download(url, fn)
            catch
                error("fail to download $fn")
            end
        end
    end            
end
