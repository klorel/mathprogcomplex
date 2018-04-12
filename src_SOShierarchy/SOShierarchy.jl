using DataStructures
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))

## Symetries
abstract type AbstractSymetry end
type PhaseInvariance <: AbstractSymetry end



mutable struct RelaxationContext
    ismultiordered
    issparse
    symmetries::OrderedSet{DataType} # ::OrderedSet{DataType}
    hierarchykind             # :Complex or :Real
    renamevars                # Replace variables with by shorter named ones
    di
    ki
end


"""
    MomentMatrix(mm, vars, order)

    Store a moment or localizing matrix of size `order`, corresponding to the `vars` variables in the `mm` dictionnary.
    **Note** that the matrix is indexed by a tuple of exponents, *the first of which contains only conjugated variables*, et second only real ones.
"""
mutable struct MomentMatrix
    mm::SortedDict{Tuple{Exponent, Exponent}, AbstractPolynomial}
    vars::OrderedSet{Variable}
    order::Int
end

"""
    momentrel = MomentRelaxationPb(obj, cstrs)

    Store a Moment Relaxation problem.
"""
struct MomentRelaxationPb
    objective::AbstractPolynomial
    constraints::SortedDict{Tuple{String, String}, MomentMatrix}
end

# """
#     MomentMatrixBasis(basis, expo2int, int2expo, msize)

#     Store the matrix coefficients of the moment variable decomposition of the moment matrix.

#     Arguments:
#     - basis::SortedDict{Exponent, AbstractMatrix} : matrice correspondant au moment clé,
#     - expo2int : carte donnant les coordonnées de l'exposant clé dans la matrice des moments initiale,
#     - int2expo : carte donnant l'xposant correspondant aux coordonnées clé dans la matrice des moments initiale,
#     - msize: taille de la matrice des moments initiale (ordre d-k).
# """
# mutable struct MomentMatrixBasis
#     basis::SortedDict{Tuple{Exponent, Exponent}, AbstractMatrix} # Les exposants sont de degré inférieur à d
#     expo2int::SortedDict{Exponent, Int}   # Carte de la matrice des coefficients d-ki
#     int2expo::SortedDict{Int, Exponent}
#     msize::Int                      # size of the matrix
# end


const SDPBody = SortedDict{Tuple{String, String, Exponent, Exponent}, SortedDict{Tuple{Exponent, Exponent}, Number}}
const SDPRhs = SortedDict{Tuple{Exponent, Exponent}, Number}

"""
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end


include("setproblem.jl")
include("momentmatrix.jl")
include("sparsity.jl")
include("symmetries.jl")
include("SDPcontainer.jl")

include("example_problems.jl")
include("utils.jl")
# include("compute_Bi.jl")
# include("build_SDP_SOS.jl")
# include("export_JuMP.jl")


function print_cmat(mat::AbstractArray, round = 1e-3)
    for i=1:size(mat, 1)
        for j=1:size(mat, 2)
            re, im = real(mat[i, j]), imag(mat[i, j])
            @printf("% 5.4f", re)
            @printf(" ")
            @printf("%+5.4fim", im)
            @printf("   ")
        end
        @printf("\n")
    end
end