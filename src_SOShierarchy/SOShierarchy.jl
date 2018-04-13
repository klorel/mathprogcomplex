using DataStructures
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))

## Symetries
abstract type AbstractSymetry end
type PhaseInvariance <: AbstractSymetry end



mutable struct RelaxationContext
    ismultiordered
    issparse
    symmetries::SortedSet{DataType} # ::SortedSet{DataType}
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
    vars::SortedSet{Variable}
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