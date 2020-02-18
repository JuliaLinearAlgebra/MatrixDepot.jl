"""
    user-defined matrix generators
    to be populated with `include_generator` in user source code
"""
const USERMATRIXDICT = Dict{String,Function}()

"""
Associate names with matrix-generating functions
"""
const MATRIXDICT = Dict("hilb" => hilb, "hadamard" => hadamard,
                  "cauchy" => cauchy, "circul" => circul,
                  "dingdong" => dingdong, "frank" => frank,
                  "invhilb" => invhilb, "forsythe" => forsythe,
                  "magic" => magic, "grcar" => grcar,
                  "triw" => triw, "moler" => moler,
                  "pascal" => pascal, "kahan" => kahan,
                  "pei" => pei, "vand" => vand,
                  "invol" => invol, "chebspec" => chebspec,
                  "lotkin" => lotkin, "clement" => clement,
                  "fiedler" => fiedler, "minij" => minij,
                  "binomial" => binomialm, "tridiag" => tridiag,
                  "lehmer" => lehmer, "parter" => parter,
                  "chow" => chow, "randcorr" => randcorr,
                  "poisson" => poisson, "neumann" => neumann,
                  "rosser" => rosser, "sampling" => sampling,
                  "wilkinson" => wilkinson, "rando" => rando,
                  "randsvd" => randsvd, "rohess" => rohess,
                  "kms" => kms, "wathen" => wathen,
                  "oscillate" => oscillate, "toeplitz" => toeplitz,
                  "hankel" => hankel, "golub" => golub,
                  "companion" => companion,
                  "prolate" => prolate, "deriv2" => deriv2,
                  "shaw" => shaw, "wing" => wing,
                  "foxgood" => foxgood, "heat" => heat,
                  "baart" => baart, "phillips" => phillips,
                  "gravity" => gravity, "blur" => blur,
                  "spikes" => spikes, "ursell" => ursell,
                  "parallax" => parallax, "erdrey" => erdrey,
                  "gilbert" => gilbert, "smallworld" => smallworld
 )

"""
    predefined matrix classes (for the generated functions)
"""
const MATRIXCLASS = Dict(
            :symmetric => ["hilb", "cauchy", "circul", "dingdong",
                             "invhilb", "moler", "pascal", "pei",
                             "clement", "fiedler", "minij", "tridiag",
                             "lehmer", "randcorr", "poisson", "wilkinson",
                             "kms", "wathen", "oscillate", "prolate",
                             "hankel"],

            :inverse => ["hilb", "hadamard", "cauchy", "invhilb",
                           "forsythe", "magic", "triw", "moler", "pascal",
                           "kahan", "pei", "vand", "invol", "lotkin",
                           "clement", "fiedler", "minij", "tridiag",
                           "lehmer", "poisson", "kms" ],

            :illcond => ["hilb", "cauchy", "frank", "invhilb",
                            "forsythe", "triw", "moler", "pascal",
                            "kahan","pei", "vand", "invol", "lotkin",
                            "tridiag", "rosser", "randsvd", "kms",
                            "oscillate", "prolate", "golub"],

            :posdef => ["hilb", "cauchy", "circul", "invhilb",
                           "moler", "pascal", "pei", "minij", "tridiag",
                           "lehmer", "poisson", "kms", "wathen", "oscillate"],

            :eigen =>   ["hadamard", "circul", "dingdong", "frank",
                           "forsythe", "grcar", "pascal", "invol","chebspec",
                           "lotkin", "clement", "fiedler", "minij",
                           "tridiag", "parter", "chow", "poisson", "neumann",
                           "rosser", "sampling", "wilkinson","wathen",
                           "oscillate"],

             :sparse => ["poisson", "neumann", "wathen", "blur", "erdrey", "gilbert",
                          "smallworld"],

             :random => ["rosser", "rando", "randcorr", "randsvd", "rohess",
                          "wathen", "oscillate", "golub", "erdrey", "gilbert", "smallworld"],

             :regprob => ["deriv2", "shaw", "wing", "foxgood", "heat",
                           "baart", "phillips", "gravity", "blur",
                           "spikes", "ursell", "parallax"],

              :graph => ["erdrey", "gilbert", "smallworld"]
)


# remote parameters for several data sources
const SS_REMOTE = SSRemoteType(RemoteParametersNew(
                    "https://sparse.tamu.edu",
                    "https://sparse.tamu.edu/MM",
                    "https://sparse.tamu.edu/files/ss_index.mat",
                    "https://sparse.tamu.edu/files/ssstats.csv",
                    ".tar.gz"
                   ))

const TA_REMOTE = TURemoteType(RemoteParameters(
                    "https://sparse.tamu.edu",
                    "https://sparse.tamu.edu/MM",
                    "https://sparse.tamu.edu/?per_page=All",
                    """<title>SuiteSparse Matrix Collection</title>""",
                    ("", """">Matrix Market""", 4, ".tar.gz", 2,
                     r"<td class='column-([[:alpha:]_]+)[^']*'>([^<]*)</td>"),
                    ".tar.gz"
                   ))

const UF_REMOTE = TURemoteType(RemoteParameters(
                    "https://www.cise.ufl.edu/research/sparse",
                    "https://www.cise.ufl.edu/research/sparse/MM",
                    "https://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html",
                    """<title>UF Sparse Matrix Collection - sorted by id</title>""",
                    ("", """>MM</a>""", 4, ".tar.gz", 2, r"<td>([[:digit:]]*)</td>"),
                    ".tar.gz"
                   ))

const MM_REMOTE = MMRemoteType(RemoteParameters(
                    "http://math.nist.gov/MatrixMarket",
                    "ftp://math.nist.gov/pub/MatrixMarket2",
                    "http://math.nist.gov/MatrixMarket/matrices.html",
                    """<TITLE>The Matrix Market Matrices by Name</TITLE>""",
                    ("M", """<A HREF="/MatrixMarket/data/""", 2, ".html", 3, nothing),
                    ".mtx.gz"
                   ))

# preferred remote source for UF matrix collection (TAMU vs. UFl)
uf_remote = SS_REMOTE # may be altered
preferred(::Type{<:Union{TURemoteType,SSRemoteType}}) = uf_remote
preferred(::Type{MMRemoteType}) = MM_REMOTE
alternate(::Type{TURemoteType}) = uf_remote === SS_REMOTE ? UF_REMOTE : SS_REMOTE
alternate(::Type{MMRemoteType}) = MM_REMOTE
toggle_remote(T=TURemoteType) = begin global uf_remote = alternate(T) end
"""
    The place to store all matrix data in process
"""
const MATRIX_DB = MatrixDatabase()

# local storage directory
const DATA_DIR = abspath(dirname(@__FILE__),"..", "data")
const MY_DEPOT_DIR = abspath(dirname(@__FILE__), "..", "myMatrixDepot")
data_dir() = get(ENV, "MATRIXDEPOT_DATA", DATA_DIR)
user_dir() = get(ENV, "MATRIXDEPOT_MYDEPOT", MY_DEPOT_DIR)
url_redirect() = URL_REDIRECT[] = get(ENV, "MATRIXDEPOT_URL_REDIRECT", "0") != "0"

const REDIRECT_DIR = abspath(dirname(@__FILE__), "..", "test", "data")
const URL_REDIRECT = Ref(false)
function redirect(url::AbstractString)
    if URL_REDIRECT[]
        urlpart = split(url, ":/", limit=2)[2]
        if Sys.iswindows()
            string("file:\\", REDIRECT_DIR, replace(urlpart, '/' => '\\'))
        else
            string("file://", REDIRECT_DIR, urlpart)
        end
    else
        url
    end
end
