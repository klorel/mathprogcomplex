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
set_vartypes!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)

Input all matrix varibales structural information regarding name and kind in the appropriate attributes of `SDP_Problem`.
"""
function set_vartypes!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  n_sdp, n_sym = 0, 0
  for i=1:size(instance.VAR_TYPES, 1)
    (block_name, block_type) = instance.VAR_TYPES[i, :]

    if block_type == "SDP"
      n_sdp += 1
      block = SDP_Block(n_sdp, block_name)
      sdp.name_to_sdpblock[block_name] = block
      sdp.id_to_sdpblock[n_sdp] = block

    elseif block_type == "Sym"
      n_sym += 1
      block = Sym_Block(n_sym, block_name)
      sdp.name_to_symblock[block_name] = block
      sdp.id_to_symblock[n_sym] = block
    else
      error("set_vartypes()!: Unknown blockvar type")
    end
  end

  if debug
    warn("TBD")
  end
end


"""
set_blocks!(sdp::SDP_Problem, instance::SDP_Instance)

Fill all variables blocks with their elementary variables previously declared by `set_vartypes!`.
"""
function set_blocks!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i in 1:size(instance.BLOCKS, 1)
    block_name, var1, var2 = instance.BLOCKS[i, 2:4]

    if haskey(sdp.name_to_sdpblock, block_name)
      cur_blockvar = sdp.name_to_sdpblock[block_name]

      # Adding vars and ids to SDP block
      if !haskey(cur_blockvar.var_to_id, var1)
        cur_blockvar.var_to_id[var1] = length(cur_blockvar.var_to_id) + 1
      end
      if !haskey(cur_blockvar.var_to_id, var2)
        cur_blockvar.var_to_id[var2] = length(cur_blockvar.var_to_id) + 1
      end
    elseif haskey(sdp.name_to_symblock, block_name)
      cur_blockvar = sdp.name_to_symblock[block_name]

      # Adding var pair to SDP block
      sortedtuple = max(var1, var2), min(var1, var2)
      if !haskey(cur_blockvar.varpairs_to_id, sortedtuple)
        cur_blockvar.varpairs_to_id[sortedtuple] = -1
      end
    else
      error("set_blocks!(): Unknown block_kind $(sdp.block_to_kind) for i=$i")
    end
  end

  # Associating id to each variable pair
  n_sym = 0
  for (blockname, block) in sdp.name_to_symblock
    for varpair in keys(block.varpairs_to_id)
      n_sym += 1
      block.varpairs_to_id[varpair] = n_sym
    end
  end

  sdp.n_scalvarsym = n_sym
  if debug
  end
end


function set_matrices!(sdp::SDP_Problem, instance::SDP_Instance; debug=false)
  for i=1:size(instance.BLOCKS, 1)
    (ctr_name, block_name, var1, var2, coeff) = instance.BLOCKS[i, :]

    # Sort variables for triangular matrix storage
    var1, var2 = min(var1, var2), max(var1, var2)

    if haskey(sdp.name_to_sdpblock, block_name)
      if !haskey(sdp.matrices, (ctr_name, block_name, var1, var2))
        sdp.matrices[(ctr_name, block_name, var1, var2)] = parse(coeff)
      else
        warn("set_matrices!(): sdp.matrices already has key ($ctr_name, $block_name, $var1, $var2) with val $(sdp.matrices[(ctr_name, block_name, var1, var2)]), $(prase(coeff))")
      end

    elseif haskey(sdp.name_to_symblock, block_name)
      if !haskey(sdp.matrices, (ctr_name, block_name, var1, var2))
        sdp.lin_matsym[(ctr_name, block_name, var1, var2)] = parse(coeff)
      else
        warn("set_matrices!(): sdp.matrices already has key ($ctr_name, $block_name, $var1, $var2) with val $(sdp.matrices[(ctr_name, block_name, var1, var2)]), $(prase(coeff))")
      end

    else
      error("set_matrices!(): Unhandled matrix var $block_name")
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

function print(io::IO, sdp::SDP_Problem)
  for (cstr, block) in sdp.name_to_sdpblock
    println(io, "  sdp   : $cstr -> $block")
  end

  for (cstr, block) in sdp.name_to_symblock
    println(io, "  sym   : $cstr -> $block")
  end

  println(io, "  objk  : $(obj_key())")

  for (name, ctr) in sdp.name_to_ctr
    println(io, "  ctr   : $name \t $ctr")
  end

  for (name, mat) in sdp.matrices
    println(io, "  matrix: $name \t $mat")
  end

  for (name, lin) in sdp.linear
    println(io, "  lin   : $name \t $lin")
  end

  for (name, cst) in sdp.cst_ctr
    println(io, "  cst   : $name \t $cst")
  end

  for (name, cst) in sdp.scalar_vars_sym
    println(io, "  sym v : $name \t $cst")
  end

  for (name, cst) in sdp.scalar_vars_ctr
    println(io, "  multv : $name \t $cst")
  end
end