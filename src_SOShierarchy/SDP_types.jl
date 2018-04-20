type SDP_Instance
  VAR_TYPES
  BLOCKS
  LINEAR
  CONST
end

function read_SDPInstance(path::String)
  BLOCKS = readdlm(joinpath(path, "blocks.sdp"), String)
  if filesize(joinpath(path, "lin.sdp")) > 0
    LINEAR = readdlm(joinpath(path, "lin.sdp"), String)
  else
    LINEAR = []
  end
  CONST = readdlm(joinpath(path, "const.sdp"), String)
  VAR_TYPES = readdlm(joinpath(path, "types.sdp"), String)

  SDP_Instance(VAR_TYPES,
               BLOCKS,
               LINEAR,
               CONST)
end

type SDP_Block
  id::Int64
  name::String
  var_to_id::Dict{String, Int64}

  SDP_Block(id::Int64, name::String) = new(id, name, Dict{String, Int64}())
  SDP_Block() = new(-1, "", Dict{String, Int64}())
end


type SDP_Problem
  name_to_block::Dict{String, SDP_Block}
  id_to_block::Dict{Int64, SDP_Block}

  name_to_ctr::Dict{String, Tuple{Int64, String, Float64, Float64}} # Id, type et bornes des contraintes
  id_to_ctr::Dict{Int64, String}

  matrices::Dict{Tuple{String, String, String, String}, Float64} # Matrices du corps des contraintes / objectif
  linear::Dict{Tuple{String, String}, Float64} # Matrice portant les parties linéaires des contraintes
  cst_ctr::Dict{String, Float64} # Constante du corp des contraintes


  SDP_Problem() = new(Dict{String, SDP_Block}(),
                      Dict{Int64, SDP_Block}(),
                      Dict{String, Tuple{Int64, String, Float64, Float64}}(),
                      Dict{Int64, String}(),
                      Dict{Tuple{String, String, String, String}, Float64}(),
                      Dict{String, Dict{String, Float64}}(),
                      Dict{String, Float64}()
  ) 
end


"""
set_constraints!(sdp::SDP_Problem, instance::SDP_Instance)

Build `name_to_ctr` with explicit constraint parameters from instance with ==0 as default.
"""
function set_constraints!(sdp::SDP_Problem, instance::SDP_Instance)
  # Collect constraints names
  ctr_names = Set{String}([instance.BLOCKS[i, 1] for i=1:size(instance.BLOCKS, 1)])
  union!(ctr_names, [instance.LINEAR[i, 1] for i=1:size(instance.LINEAR, 1)])
  union!(ctr_names, [instance.CONST[i, 1] for i=1:size(instance.CONST, 1)])

  @show ctr_names

  @show sdp.name_to_ctr

  # Default contraint is EQ, == 0 
  # TODO: is that ok ?
  # TODO : should there be sparse storage here ?
  ctr_id = 1
  for ctr_name in ctr_names
    sdp.name_to_ctr[ctr_name] = (ctr_id, "EQ", 0, 0)
    sdp.id_to_ctr[ctr_id] = ctr_name
    ctr_id += 1
  end

  # # Looping over explicitly defined constraints
  # for i=1:size(instance.CONSTRAINTS, 1)
  #   (ctr_name, rhs) = instance.CONSTRAINTS[i, :]
    
  #   id_ctr = sdp.name_to_ctr[ctr_name][1]
  #   sdp.name_to_ctr[ctr_name] = (id_ctr, "EQ", parse(rhs), parse(rhs))
  # end
end


"""
set_blocks!(sdp::SDP_Problem, instance::SDP_Instance)

Set `name_to_block` with all variable blocks.
"""
function set_blocks!(sdp::SDP_Problem, instance::SDP_Instance)
  for i in 1:size(instance.BLOCKS, 1)
    block_name, var1, var2 = instance.BLOCKS[i, 2:4]
    
    if !haskey(sdp.name_to_block, block_name)
      block_id = 1 + length(sdp.name_to_block)
      sdp_block = SDP_Block(block_id, block_name)
      sdp.name_to_block[block_name] = sdp_block
      sdp.id_to_block[block_id] = sdp_block
    else
      sdp_block = sdp.name_to_block[block_name]
    end

    # Adding vars to SDP block
    if !haskey(sdp_block.var_to_id, var1)
      sdp_block.var_to_id[var1] = length(sdp_block.var_to_id) + 1
    end
    if !haskey(sdp_block.var_to_id, var2)
      sdp_block.var_to_id[var2] = length(sdp_block.var_to_id) + 1
    end
  end
end


function set_matrices!(sdp::SDP_Problem, instance::SDP_Instance)
  for i=1:size(instance.BLOCKS, 1)
    (ctr_name, block_name, var1, var2, coeff) = instance.BLOCKS[i, :]

    @assert !haskey(sdp.matrices, (ctr_name, block_name, var1, var2))

    # Sort variables for triangular matrix storage
    var1, var2 = min(var1, var2), max(var1, var2)
    sdp.matrices[(ctr_name, block_name, var1, var2)] = parse(coeff)
  end
end

function set_linear!(sdp::SDP_Problem, instance::SDP_Instance)
  # TODO
  if length(instance.LINEAR) != 0
    warn("set_linmatrix!() not implemented yet. Passing.")
  end
end

function set_const!(sdp::SDP_Problem, instance::SDP_Instance)
  for i=1:size(instance.CONST, 1)
    (ctr_name, coeff) = instance.CONST[i, :]

    @assert !haskey(sdp.cst_ctr, ctr_name)
    sdp.cst_ctr[ctr_name] = parse(coeff)
  end
end