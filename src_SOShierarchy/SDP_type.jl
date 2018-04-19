
type SDP_Instance
  ##############################################################################
  name::String
  formulation::String
  mean::String
  ##############################################################################
  BLOCKS::Matrix{String}
  CONSTRAINTS::Matrix{Any}
  MATRICES::Matrix{Any}
  CONSTANTS::Matrix{Any}
  ##############################################################################
  blocks
  constraints
  matrices
  constants::Matrix{Any}
  constants_per_constraint::Dict{String,Float64}
  DIFFERENT_BLOCKS::Set{String}
  block_var::Dict{String, Set{String}}
  ##############################################################################
  # jumpX::Array{Matrix{JuMP.Variable}}
  coeff_block::Dict{Tuple{String,String}, Set{String}}
  # mat_var::Dict{Tuple{String,String}, Dict{String,JuMP.Variable}}
  ##############################################################################
  SDP_Instance() = new()
  ##############################################################################
end

type SDP_Block
  id::Int64
  name::String
  var_to_id::Dict{String, Int64}
  SDP_Block(id::Int64) =
    new(
      Int64(id),
      String(""),
      Dict{String, Int64}()
    )
  SDP_Block() =
      new(
        Int64(-1),
        String(""),
        Dict{String, Int64}()
      )
end

type SDP_Problem
  name_to_block::Dict{String, SDP_Block}
  id_to_block::Dict{Int64, SDP_Block}

  name_to_kind::Dict{String, String}
  name_to_ctr::Dict{String, Tuple{Int64, String, Float64, Float64}}
  id_to_ctr::Dict{Int64, String}

  term_to_block::Dict{Tuple{String, String}, SDP_Block}
  term_to_all_block::Dict{Tuple{String, String}, Set{String}}
  matrix::Dict{String, Dict{Tuple{String,String,String}, Float64}}

  cst_ctr::Dict{String, Float64}

  SDP_Problem() =
    new(
      Dict{String, SDP_Block}(),
      Dict{Int64, SDP_Block}(),
      Dict{String, String}(),
      Dict{String, Tuple{Int64, String, Float64, Float64}}(),
      Dict{Int64, String}(),
      Dict{Tuple{String, String}, SDP_Block}(),
      Dict{Tuple{String, String}, Dict{String, SDP_Block}}(),
      Dict{String, Dict{Tuple{String,String,String}, Float64}}(),
      Dict{String, Float64}(),
      )
end

function add_blocks(instance::SDP_Instance, sdp::SDP_Problem)
  for i in 1:size(instance.BLOCKS)[1]
    # println(instance.BLOCKS[i, :])
    block = instance.BLOCKS[i, 1]
    var =  instance.BLOCKS[i, 2]
    # println("var   = ", var)
    # println("block = ", block)
    if !haskey(sdp.name_to_block, block)
      block_id = Int64(1 + length(sdp.name_to_block))
      sdp_block = SDP_Block(block_id)
      sdp_block.name = block
      sdp.name_to_block[block] = sdp_block
      sdp.id_to_block[block_id] = sdp_block
    else
      sdp_block = sdp.name_to_block[block]
    end
    if !haskey(sdp_block.var_to_id, var)
      var_id = Int64(1+length(sdp_block.var_to_id))
      sdp_block.var_to_id[var] =  var_id
    end
  end
  for name_block in sdp.name_to_block
    name = name_block[1]
    sdp_block = name_block[2]
    sdp.matrix[name] = Dict{Tuple{String,String}, Float64}()
    for var1 in keys(sdp_block.var_to_id)
      for var2 in keys(sdp_block.var_to_id)
        if !haskey(sdp.term_to_block, (var1, var2))
          sdp.term_to_block[var1, var2] = sdp_block
          sdp.term_to_all_block[var1, var2] = Set{String}([name])
        else
          push!(sdp.term_to_all_block[var1, var2], name)
          if length(sdp.term_to_block[var1, var2].var_to_id)>length(sdp_block.var_to_id)
            sdp.term_to_block[var1, var2] = sdp_block
          end
        end
        # println(var1, " - ", var2, " : ", sdp.term_to_all_block[var1, var2])
      end
    end
    # println(name)
    # println(keys(sdp_block.var_to_id))
  end
end

function add_constraints(instance::SDP_Instance, sdp::SDP_Problem)
  # lowerbound and upperbound
  for i in 1:size(instance.CONSTRAINTS)[1]
    (name, kind, lb, ub)  = instance.CONSTRAINTS[i, 1:4]
    # kind = name_kind_lb_ub[2]
    # lb = name_kind_lb_ub[3]
    # ub  = name_kind_lb_ub[4]
    sdp.name_to_kind[name] = kind
    # @printf("%20s%20s%20s%20s\n", name, kind, lb, ub)

    ctr_id = 1+length(sdp.name_to_ctr)
    sdp.name_to_ctr[name] = (ctr_id, kind, lb*C_SCALING_FACTOR(), ub*C_SCALING_FACTOR())
    sdp.cst_ctr[name] = 0

    # if kind != "RNG"
    #   ctr_id = 1+length(sdp.name_to_ctr)
    #   sdp.name_to_ctr[name] = (ctr_id, kind, lb)
    #   sdp.cst_ctr[name] = 0
    # else
    #   @assert(kind == "RNG")
    #
    #   ctr_id = 1+length(sdp.name_to_ctr)
    #   sdp.name_to_ctr[name*"_GEQ"] =  (ctr_id, "GEQ", lb)
    #   sdp.cst_ctr[name*"_GEQ"] = 0
    #
    #   ctr_id = 1+length(sdp.name_to_ctr)
    #   sdp.name_to_ctr[name*"_LEQ"] = (ctr_id, "LEQ", ub)
    #   sdp.cst_ctr[name*"_LEQ"] = 0
    # end
  end
  sdp.cst_ctr["OBJ"] = 0
  # for k in keys(sdp.cst_ctr)
  #   println(k)
  # end
end

function add_matrix(instance::SDP_Instance, sdp::SDP_Problem)
  for i in 1:size(instance.MATRICES)[1]
    (name, var1, var2, coeff)  = instance.MATRICES[i, 1:4]
    if var1!="NONE" && var2!="NONE"
      block_name = sdp.term_to_block[var1, var2].name

      if name != "OBJ"
        kind = sdp.name_to_kind[name]

        if !haskey(sdp.matrix, name)
          sdp.matrix[name] = Dict{Tuple{String,String,String}, Float64}()
        end
        sdp.matrix[name][block_name, min(var1, var2), max(var1, var2)] = coeff*X_SCALING_FACTOR()*X_SCALING_FACTOR()*C_SCALING_FACTOR()


        # if kind != "RNG"
        #   if !haskey(sdp.matrix, name)
        #     sdp.matrix[name] = Dict{Tuple{String,String,String}, Float64}()
        #   end
        #   sdp.matrix[name][block_name, min(var1, var2), max(var1, var2)] = coeff
        # else
        #   if !haskey(sdp.matrix, name*"_GEQ")
        #     sdp.matrix[name*"_GEQ"] = Dict{Tuple{String,String,String}, Float64}()
        #     sdp.matrix[name*"_LEQ"] = Dict{Tuple{String,String,String}, Float64}()
        #   end
        #   sdp.matrix[name*"_GEQ"][block_name, min(var1, var2), max(var1, var2)] = coeff
        #   sdp.matrix[name*"_LEQ"][block_name, min(var1, var2), max(var1, var2)] = coeff
        # end
      else
        if !haskey(sdp.matrix, name)
          sdp.matrix[name] = Dict{Tuple{String,String,String}, Float64}()
        end
        sdp.matrix[name][block_name, min(var1, var2), max(var1, var2)] = coeff*X_SCALING_FACTOR()*X_SCALING_FACTOR()*C_SCALING_FACTOR()
      end
    else
      @assert(var1=="NONE" && var2=="NONE")
      sdp.cst_ctr[name] += coeff*C_SCALING_FACTOR()
    end
  end
end

function add_coupling(instance::SDP_Instance, sdp::SDP_Problem)
  for term_all_block in sdp.term_to_all_block
    var1 = term_all_block[1][1]
    var2 = term_all_block[1][2]
    all_block = [s for s in term_all_block[2]]
    # println(all_block)
    i_chosen = 1
    min_size = length(sdp.name_to_block[all_block[i_chosen]].var_to_id)
    for i in 1:length(all_block)
      if all_block[i] == sdp.term_to_block[var1,var2].name
        i_chosen = i
      end
    end
    for i in 1:length(all_block)
      if i != i_chosen
        block_1 = all_block[i_chosen]
        block_2 = all_block[i]
        ctr_name = "cc_"*var1*"_"*var2*"_"*block_1*"_"*block_2

        ctr_id = 1+length(sdp.name_to_ctr)
        sdp.name_to_ctr[ctr_name] = (ctr_id, "EQ", 0, 0)
        sdp.matrix[ctr_name] = Dict{Tuple{String,String,String}, Float64}()
        sdp.matrix[ctr_name][block_1, min(var1, var2), max(var1, var2)] = +1
        sdp.matrix[ctr_name][block_2, min(var1, var2), max(var1, var2)] = -1
        sdp.cst_ctr[ctr_name] = 0
        # @printf("term is %15s %15s coupling %6s and %6s in %15s\n", var1, var2, block_1, block_2, ctr_name)
      end
    end
  end
end

function init(instance::SDP_Instance)
  instance.name = ""
  instance.formulation = "DENSE"
  instance.mean = "NO"
  ##############################################################################
  instance.BLOCKS = Matrix{String}(0, 0)
  instance.CONSTRAINTS = Matrix{Any}(0, 0)
  instance.MATRICES = Matrix{Any}(0, 0)
  instance.CONSTANTS = Matrix{Any}(0, 0)
  # instance.MATRICES_without_CONSTANTS = Matrix{Any}()
  ##############################################################################
  instance.constants = Matrix{Any}(0, 0)
  instance.constants_per_constraint = Dict{String,Float64}()
  instance.DIFFERENT_BLOCKS = Set{String}()
  instance.block_var = Dict{String, Set{String}}()
  ##############################################################################
  # instance.jumpX=Array{Matrix{JuMP.Variable}}(NB_BLOCKS)
  # instance.coeff_block = Dict{Tuple{String,String}, Set{String}}()
  # instance.mat_var = Dict{Tuple{String,String}, Dict{String,JuMP.Variable}}()
  ##############################################################################
end
const MatrixSDP = Dict{Tuple{String,String}, Float64}
const Column = Dict{String, MatrixSDP};

type ColumnGeneration
  problem::SDP_Problem
end
