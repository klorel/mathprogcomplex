using DataStructures
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))

## Symetries
abstract type AbstractSymetry end
type PhaseInvariance <: AbstractSymetry end



mutable struct RelaxationContext
    ismultiordered
    issparse
    symmetries::SortedSet{DataType} # ::SortedSet{DataType}
    hierarchykind                   # :Complex or :Real
    renamevars                      # Replace variables with by shorter named ones
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
    matrixkind::Symbol            # Either :SDP or :Sym
end

"""
    momentrel = MomentRelaxationPb(obj, cstrs)

    Store a Moment Relaxation problem.
"""
struct MomentRelaxationPb
    objective::AbstractPolynomial
    constraints::SortedDict{Tuple{String, String}, MomentMatrix}
    vars_overlap::SortedDict{Variable, SortedSet{String}}
end



# type MatVar
#   name::String
#   kind::Symbol # Either :SDP or :Sym
# end

# mutable struct SDPInstance_
#     sdpkind::Symbol # Either :Complex or :Real

#     name_to_sdpblock::Dict{String, MatVar}
#     id_to_sdpblock::Dict{Int64, MatVar}

#     name_to_ctr::SortedDict{Tuple{Exponent, Exponent}, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
#     id_to_ctr::SortedDict{Int64, Tuple{Exponent, Exponent}}
#     obj_name::String

#     matrices::SortedDict{Tuple{Tuple{Exponent, Exponent}, String, Exponent, Exponent}, Number}
#     # ((α, β), block_name, γ, δ) -> coeff
#     lin::SortedDict{Tuple{Tuple{Exponent, Exponent}, Exponent}, Number}
#     # ((α, β), var) -> coeff
#     cst_ctr::SortedDict{Tuple{Exponent, Exponent}, Number}
#     # (α, β) -> coeff

#     SDPInstance() = new(:None,
#                         Dict{String, MatVar}(),
#                         Dict{Int64, MatVar}(),
#                         SortedDict{Tuple{Exponent, Exponent}, Tuple{Int64, String, Float64, Float64}}(), # Id, type et bornes des contraintes
#                         SortedDict{Int64, Tuple{Exponent, Exponent}}(),
#                         "",
#                         SortedDict{Tuple{Tuple{Exponent, Exponent}, String, Exponent, Exponent}, Number}(),
#                         SortedDict{Tuple{Tuple{Exponent, Exponent}, Exponent}, Number}(),
#                         SortedDict{Tuple{Exponent, Exponent}, Number}()
#                         )
# end

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

type Var_Block
  id::Int64
  name::String
  var_to_id::SortedDict{String, Int64}
  varkind::Symbol

  Var_Block(id::Int64, name::String, varkind) = new(id, name, SortedDict{String, Int64}(), varkind)
end


type SDP_Problem
  # Matrix variables
  block_to_kind::SortedDict{String, Symbol}

  # SDP vars
  name_to_sdpblock::SortedDict{String, Var_Block}
  id_to_sdpblock::SortedDict{Int64, Var_Block}

  # Free vars
  name_to_symblock::SortedDict{String, Var_Block}
  id_to_symblock::SortedDict{Int64, Var_Block}

  name_to_ctr::SortedDict{String, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
  id_to_ctr::SortedDict{Int64, String}
  obj_name::String

  matrices::SortedDict{Tuple{String, String, String, String}, Float64} # Matrices du corps des contraintes / objectif
  linear::SortedDict{Tuple{String, String, String, String}, Float64} # Matrice portant les parties linéaires des contraintes
  cst_ctr::SortedDict{String, Float64} # Constante du corp des contraintes

  scalar_vars_sym::Dict{Tuple{String, String, String}, Int64}
  scalar_vars_ctr::Dict{Tuple{String, String}, Int64}

  SDP_Problem() = new(SortedDict{String, Symbol}(),
                      SortedDict{String, Var_Block}(),
                      SortedDict{Int64, Var_Block}(),
                      SortedDict{String, Var_Block}(),
                      SortedDict{Int64, Var_Block}(),
                      SortedDict{String, Tuple{Int64, String, Float64, Float64}}(),
                      SortedDict{Int64, String}(),
                      obj_key(),
                      SortedDict{Tuple{String, String, String, String}, Float64}(),
                      SortedDict{Tuple{String, String, String, String}, Float64}(),
                      SortedDict{String, Float64}(),
                      Dict{Tuple{String, String, String}, Int64}(),
                      Dict{Tuple{String, String}, Int64}()
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
include("run_hierarchy.jl")

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