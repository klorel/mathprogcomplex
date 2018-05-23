using DataStructures
include(joinpath(ROOT, "src_PolynomialOptim", "PolynomialOptim.jl"))


###############################################################################
## Relaxation context, symmetries and cliques
###############################################################################
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
    SparsityPattern

    Type for storing and working on sparsitty patterns.
"""
type SparsityPattern end

include("build_relctx.jl")
include("build_maxcliques.jl")

abstract type AbstractSymetry end
type PhaseInvariance <: AbstractSymetry end
include("symmetries.jl")

###############################################################################
## Moment Problem
###############################################################################
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
    vars_overlap::SortedDict{Exponent, SortedSet{String}}
end

include("build_momentpb.jl")



###############################################################################
## SOS Problem
###############################################################################
const SDPBlocks = SortedDict{Tuple{Tuple{Exponent, Exponent}, String, Exponent, Exponent}, Number} # ((α, β), block_name, γ, δ) -> coeff
const SDPLinSym = SortedDict{Tuple{Tuple{Exponent, Exponent}, String, Exponent}, Number}           # ((α, β), block_name, var) -> coeff
const SDPLin = SortedDict{Tuple{Tuple{Exponent, Exponent}, Exponent}, Number}                      # ((α, β), var) -> coeff
const SDPCst = SortedDict{Tuple{Exponent, Exponent}, Number}                                       # (α, β) -> coeff

mutable struct SDPInstance
    block_to_vartype::SortedDict{String, Symbol}  # Either :SDP, :Sym
    blocks::SDPBlocks
    linsym::SDPLinSym
    lin::SDPLin
    cst::SDPCst
end

include("build_SDPInstance.jl")
include("export_SDPInstance.jl")
include("SDPInstance_cplx2real.jl")



###############################################################################
## Mosek Structures
###############################################################################
type SDP_Instance
  VAR_TYPES
  BLOCKS
  LINEAR
  CONST
end


type SDP_Block
  id::Int64
  name::String
  var_to_id::SortedDict{String, Int64}

  SDP_Block(id::Int64, name::String) = new(id, name, SortedDict{String, Int64}())
end

type SDP_Problem
  # SDP vars
  name_to_sdpblock::SortedDict{String, SDP_Block}
  id_to_sdpblock::SortedDict{Int64, SDP_Block}

  # Scalar variables
  scalvar_to_id::Dict{String, Int64}

  # Objective / constraints
  obj_name::Tuple{String, String}
  name_to_ctr::SortedDict{Tuple{String, String}, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
  id_to_ctr::SortedDict{Int64, Tuple{String, String}}

  matrices::SortedDict{Tuple{Tuple{String, String}, String, String, String}, Float64} # Matrices SDP du corps des contraintes / objectif
  linear::SortedDict{Tuple{Tuple{String, String}, String}, Float64} # Matrice portant les parties linéaires des contraintes
  cst_ctr::SortedDict{Tuple{String, String}, Float64} # Constante du corps des contraintes

  SDP_Problem() = new(SortedDict{String, SDP_Block}(),
                      SortedDict{Int64, SDP_Block}(),
                      Dict{String, Int64}(),
                      obj_key(),
                      SortedDict{Tuple{String, String}, Tuple{Int64, String, Float64, Float64}}(),
                      SortedDict{Int64, Tuple{String, String}}(),
                      SortedDict{Tuple{Tuple{String, String}, String, String, String}, Float64}(),
                      SortedDict{Tuple{Tuple{String, String}, String}, Float64}(),
                      SortedDict{Tuple{String, String}, Float64}()
                      )
end

include("run_mosek.jl")
include("build_SDP_Problem.jl")



###############################################################################
## Unsorted
###############################################################################


include("run_hierarchy.jl")

include("example_problems.jl")
include("utils.jl")