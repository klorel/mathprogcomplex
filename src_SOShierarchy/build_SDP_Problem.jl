function read_SDPInstance(path::String)
  BLOCKS = readdlm(joinpath(path, "blocks.sdp"), String)
  if isfile(joinpath(path, "lin.sdp"))
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



"""
set_constraints!(sdp::SDP_Problem, instance::SDP_Instance)

Build `name_to_ctr` with explicit constraint parameters from instance with ==0 as default.
"""
function set_constraints!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  # Collect constraints names
  ctr_names = SortedSet{String}([instance.BLOCKS[i, 1] for i=1:size(instance.BLOCKS, 1)])
  union!(ctr_names, [instance.LINEAR[i, 1] for i=1:size(instance.LINEAR, 1)])
  union!(ctr_names, [instance.CONST[i, 1] for i=1:size(instance.CONST, 1)])

  if haskey(ctr_names, obj_key())
    delete!(ctr_names, obj_key())
  else
    warn("No ctrkey matching objective key $(obj_key())")
  end

  # Default contraint is EQ, == 0
  # TODO: is that ok ?
  # TODO : should there be sparse storage here ?
  ctr_id = 1
  for ctr_name in ctr_names
    sdp.name_to_ctr[ctr_name] = (ctr_id, "EQ", 0, 0)
    sdp.id_to_ctr[ctr_id] = ctr_name
    ctr_id += 1
  end

  if debug

  end
end


"""
set_blocks!(sdp::SDP_Problem, instance::SDP_Instance)

Set `name_to_block` with all variable blocks.
"""
function set_blocks!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
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

  if debug
  end
end


function set_matrices!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.BLOCKS, 1)
    (ctr_name, block_name, var1, var2, coeff) = instance.BLOCKS[i, :]

    # Sort variables for triangular matrix storage
    var1, var2 = min(var1, var2), max(var1, var2)
    if !haskey(sdp.matrices, (ctr_name, block_name, var1, var2))
      sdp.matrices[(ctr_name, block_name, var1, var2)] = parse(coeff)
    else
      warn("set_matrices!() : sdp.matrices already has key ($ctr_name, $block_name, $var1, $var2) with val $(sdp.matrices[(ctr_name, block_name, var1, var2)]), $(prase(coeff))")
    end
  end

  if debug
  end
end

function set_linear!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  # TODO
  if length(instance.LINEAR) != 0
    warn("  set_linmatrix!() not implemented yet. Passing.")
  end

  if debug
  end
end

function set_const!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.CONST, 1)
    (ctr_name, coeff) = instance.CONST[i, :]

    @assert !haskey(sdp.cst_ctr, ctr_name)
    sdp.cst_ctr[ctr_name] = parse(coeff)
  end

  if debug
  end
end

function set_vartypes!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.VAR_TYPES, 1)
    (ctrname, vartype) = instance.VAR_TYPES[i, :]
    if vartype == "SYM"
      warn("Not enforcing symmetry of matrix $ctrname")
    end
  end
  if debug
    warn("TBD")
  end
end

function print(io::IO, sdp::SDP_Problem)
  for (cstr, block) in sdp.name_to_block
    println(io, "  b $cstr -> $block")
  end

  println(io, "  objkey $(obj_key())")

  for (name, ctr) in sdp.name_to_ctr
    println(io, "  * $name \t $ctr")
  end

  for (name, mat) in sdp.matrices
    println(io, "  s $name \t $mat")
  end

  for (name, lin) in sdp.linear
    println(io, "  l $name \t $lin")
  end

  for (name, cst) in sdp.cst_ctr
    println(io, "  c $name \t $cst")
  end
end