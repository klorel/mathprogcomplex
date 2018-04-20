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
    cstrtypes
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


# mutable struct SDPInstance
#     name
#     sense::Symbol       # :min or :max
#     kind::Symbol        # :Real or :Complex

#     objective::SDPForm
#     constraints::SortedDict{String, SDPForm}

#     isexportready::Bool
#     isconsistent::Bool

#     SDPInstance() = new("", :Undef, :Undef, 
#                         SDPForm(), SortedDict{String, SDPForm}(),
#                         false, false)
# end

# mutable struct SDPForm
#     name::String,
#     blocks::SortedDict{String, SDPBlock}
#     lin::Polynomial
#     cst::Number
#     ub::Number
#     lb::Number

#     SDPForm() = new("", SortedDict{String, SDPBlock}(), Polynomial(), NaN, NaN, NaN)
# end

const SDPBlock = SortedDict{Tuple{Exponent, Exponent}, Number}

const SDPBlocks = SortedDict{Tuple{Tuple{Exponent, Exponent}, String, Exponent, Exponent}, Number}
# ((α, β), block_name, γ, δ) -> coeff

const SDPLin = SortedDict{Tuple{Tuple{Exponent, Exponent}, Exponent}, Number}
# ((α, β), var) -> coeff

const SDPcst = SortedDict{Tuple{Exponent, Exponent}, Number}
# (α, β) -> coeff

mutable struct SDPInstance
    blocks::SDPBlocks
    lin::SDPLin
    cst::SDPcst
end



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
include("export_sdp.jl")
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