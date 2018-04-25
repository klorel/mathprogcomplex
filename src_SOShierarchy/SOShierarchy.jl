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


type SDP_Instance
  VAR_TYPES
  BLOCKS
  LINEAR
  CONST
end

type SDP_Block
  id::Int64
  name::String
  var_to_id::Dict{String, Int64}

  SDP_Block(id::Int64, name::String) = new(id, name, Dict{String, Int64}())
  SDP_Block() = new(-1, "", Dict{String, Int64}())
end


type SDP_Problem
  name_to_block::SortedDict{String, SDP_Block}
  id_to_block::SortedDict{Int64, SDP_Block}

  name_to_ctr::SortedDict{String, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
  id_to_ctr::SortedDict{Int64, String}

  matrices::SortedDict{Tuple{String, String, String, String}, Float64} # Matrices du corps des contraintes / objectif
  linear::SortedDict{Tuple{String, String}, Float64} # Matrice portant les parties linéaires des contraintes
  cst_ctr::SortedDict{String, Float64} # Constante du corp des contraintes


  SDP_Problem() = new(SortedDict{String, SDP_Block}(),
                      SortedDict{Int64, SDP_Block}(),
                      SortedDict{String, Tuple{Int64, String, Float64, Float64}}(),
                      SortedDict{Int64, String}(),
                      SortedDict{Tuple{String, String, String, String}, Float64}(),
                      SortedDict{Tuple{String, String}, Float64}(),
                      SortedDict{String, Float64}()
  )
end

"""
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end

include("build_relctx.jl")
include("build_maxcliques.jl")
include("build_momentpb.jl")
include("build_SDPInstance.jl")
include("build_SDP_Problem.jl")

include("symmetries.jl")
include("export_SDPInstance.jl")

include("example_problems.jl")
include("run_mosek.jl")
include("utils.jl")


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