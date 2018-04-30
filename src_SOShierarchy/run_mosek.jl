using Mosek

printstream(msg::String) = print(msg)

function get_triplets(problem::SDP_Problem; debug = false)

  nzc = 0
  nza = 0
  for ((objctr, block, var1, var2), coeff) in problem.matrices

      if objctr == obj_key()
        nzc += 1
      else
        nza += 1
      end
  end
  println("nza : ", nza)
  println("nzc : ", nzc)
  barai = Int32[0 for i in 1:nza]
  baraj = Int32[0 for i in 1:nza]
  barak = Int32[0 for i in 1:nza]
  baral = Int32[0 for i in 1:nza]
  baraijkl = Float64[0 for i in 1:nza]

  barcj = Int32[0 for i in 1:nzc]
  barck = Int32[0 for i in 1:nzc]
  barcl = Int32[0 for i in 1:nzc]
  barcjkl = Float64[0 for i in 1:nzc]

  nzc=0
  nza=0
  for ((objctr, block, var1, var2), coeff) in problem.matrices
    sdp_block = problem.name_to_block[block]
    lower = min(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
    upper = max(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
    if objctr == problem.obj_name
      nzc+=1
      barcj[nzc] = sdp_block.id
      barck[nzc] = upper
      barcl[nzc] = lower
      barcjkl[nzc] = coeff
    else
      nza+=1
      barai[nza] = problem.name_to_ctr[objctr][1]
      baraj[nza] = sdp_block.id
      barak[nza] = upper
      baral[nza] = lower
      baraijkl[nza] = coeff
    end
  end

  if debug
    println("*********************************************************************************")
    println("Debug -> Reading SDP_Problem")
    @printf("%5s  %5s  %5s  %s\n", "barcj", "barck", "barcl", "barcjkl")
    for i=1:length(barcj)
        @printf("%5i  %5i  %5i  %f\n", barcj[i], barck[i], barcl[i], barcjkl[i])
    end

    @printf("%5s  %5s  %5s  %5s  %s\n", "barai", "baraj", "barak", "baral", "baraijkl")
    for i=1:length(barai)
        @printf("%5i  %5i  %5i  %5i  %f\n", barai[i], baraj[i], barak[i], baral[i], baraijkl[i])
    end
  end
  return barcj, barck, barcl, barcjkl, barai, baraj, barak, baral, baraijkl
end

function get_bounds(problem::SDP_Problem; debug = false)
  numcon=length(problem.name_to_ctr)
  MOSEK_KIND = Dict(["EQ"=>MSK_BK_FX, "GEQ"=>MSK_BK_LO, "LEQ"=>MSK_BK_UP, "RNG"=>MSK_BK_RA])
  bkc = Boundkey[ Mosek.Boundkey(1)  for kv in problem.name_to_ctr]
  # constraint
  buc = Float64[0 for i in 1:numcon]
  blc = Float64[0 for i in 1:numcon]
  for (ctrname, ctr) in problem.name_to_ctr
    # @printf("%10d%20s%10d%20s\n", ctr[1], ctrname, MOSEK_KIND[ctr[2]], ctr[2])
    id_ctr=ctr[1]
    lb = ctr[3]
    ub = ctr[4]
    cst = 0
    if haskey(problem.cst_ctr, ctrname)
      cst = problem.cst_ctr[ctrname]
    end
    bkc[id_ctr] = MOSEK_KIND[ctr[2]]
    if bkc[id_ctr] == MSK_BK_UP
      buc[id_ctr] = ub - cst
    elseif bkc[id_ctr] == MSK_BK_LO
      blc[id_ctr] = lb - cst
    elseif bkc[id_ctr] == MSK_BK_FX
      blc[id_ctr] = lb - cst
      buc[id_ctr] = ub - cst
    elseif bkc[id_ctr] == MSK_BK_RA
      blc[id_ctr] = lb - cst
      buc[id_ctr] = ub - cst
    else
      error("get_bounds() : Unknown constraint kind $(ctr[2]) $(bkc[id_ctr]) $(MSK_BK_FX[1])")
    end
  end
  if debug
    warn("get_bounds(): done")
    @show numcon
    @show bkc
    @show blc
    @show buc
  end
  return numcon, bkc, blc, buc
end


function solve_mosek(problem::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64},
                                       dual::SortedDict{Tuple{String, String, String}, Float64};
                                       debug = false,
                                       logname = "")
  empty!(primal)
  empty!(dual)
  primobj = NaN
  dualobj = NaN

  num_block = length(problem.id_to_block)
  # println("num_block : ", num_block)
  barvardim = [ length(problem.id_to_block[block].var_to_id) for block in 1:num_block ]
  vardim = 0
  # println(barvardim)
  println("num_block = ",   num_block)
  numcon, bkc, blc, buc = get_bounds(problem)
  println("numcon = ",   numcon)

  barcj, barck, barcl, barcjkl , barai, baraj, barak, baral, baraijkl = get_triplets(problem, debug=debug)

  # Create a task object and attach log stream printer
  maketask() do task
      if logname != ""
        linkfiletostream(task, MSK_STREAM_LOG, logname, 0)
      else
        putstreamfunc(task,MSK_STREAM_LOG,printstream)
      end

      # Append matrix variables of sizes in 'BARVARDIM'.
      # The variables will initially be fixed at zero.
      appendbarvars(task,barvardim)
      appendvars(task, vardim)
      # Append 'numcon' empty constraints.
      # The constraints will initially have no bounds.
      appendcons(task,numcon)
      # Set the bounds on constraints.
      putconboundslice(task,1,numcon+1, bkc,blc,buc)

      # Minimize
      putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

      # Set constraints matrices
      putbarablocktriplet(task, length(barai), barai, baraj, barak, baral, baraijkl)

      # Set contraints linear part
      # putaijlist(task, ai, aj, aij)

      # Objective matrices and constant
      putbarcblocktriplet(task, length(barcj), barcj, barck, barcl, barcjkl)
      if haskey(problem.cst_ctr, problem.obj_name)
        putcfix(task, problem.cst_ctr[problem.obj_name])
      end

      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE)
      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_AGGRESSIVE)

      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_MU_RED, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-12)
      # putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-1)


      if debug
        println("*********************************************************************************")
        println("Debug -> Reading Mosek params")
        println("MSK_DPAR_INTPNT_CO_TOL_DFEAS,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_DFEAS))")
        println("MSK_DPAR_INTPNT_CO_TOL_INFEAS,   $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_INFEAS))")
        println("MSK_DPAR_INTPNT_CO_TOL_MU_RED,   $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_MU_RED))")
        println("MSK_DPAR_INTPNT_CO_TOL_PFEAS,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_PFEAS))")
        println("MSK_DPAR_INTPNT_CO_TOL_REL_GAP,  $(getdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP))")

        println("Debug -> Reading Mosek problem")
        num, subcj, subck, subcl, valcjkl = getbarcblocktriplet(task)
        @printf("%5s  %5s  %5s  %s\n", "subcj", "subck", "subcl", "valcjkl")
        for i=1:num
            @printf("%5i  %5i  %5i  %f\n", subcj[i], subck[i], subcl[i], valcjkl[i])
        end

        num, subai, subaj, subak, subal, valajkl = getbarablocktriplet(task)
        @printf("%5s  %5s  %5s  %5s  %s\n", "subai", "subaj", "subak", "subal", "valajkl")
        for i=1:num
            @printf("%5i  %5i  %5i  %5i  %f\n", subai[i], subaj[i], subak[i], subal[i], valajkl[i])
        end

        boundkeys, lbs, ubs = getconboundslice(task, 1, numcon+1)
        @show boundkeys
        @show lbs
        @show ubs
        println("*********************************************************************************")
      end

      optimize(task)
      solutionsummary(task,MSK_STREAM_MSG)

      # Get status information about the solution
      prosta = getprosta(task,MSK_SOL_ITR)
      solsta = getsolsta(task,MSK_SOL_ITR)
      if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
        # Get objective
        primobj = getprimalobj(task, MSK_SOL_ITR)
        dualobj = getdualobj(task, MSK_SOL_ITR)

        # Get primal solution
        for (id, block) in problem.id_to_block
          barx = getbarxj(task, MSK_SOL_ITR, id)
          all_variables = ["" for kv in block.var_to_id]
          for (var, varid) in block.var_to_id
            all_variables[varid] = var
          end
          n = 0
          for j in 1:length(all_variables)
            for i in j:length(all_variables)
              n+=1
              # @printf("%15s%15s%15s%25.10f\n", id_block[2].name, all_variables[i], all_variables[j], barx[n])
              primal[block.name, all_variables[i], all_variables[j]] = barx[n]
            end
          end
        end

        # Get dual solution
        for (blockid, block) in problem.id_to_block
          bars = getbarsj(task, MSK_SOL_ITR, blockid)
          id_to_var = OrderedDict([id=>var for (var,id) in block.var_to_id])
          it = 0
          for j=1:length(id_to_var), i=j:length(id_to_var)
            it += 1
            # println("($(block.name), $(id_to_var[i]), $(id_to_var[j])) = $(bars[it])")
            dual[(block.name, id_to_var[i], id_to_var[j])] = bars[it]
          end
        end

      elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
          println("Primal or dual infeasibility.\n")
      elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
          println("Primal or dual infeasibility.\n")
      elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
          println("Primal or dual infeasibility.\n")
      elseif  solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
          println("Primal or dual infeasibility.\n")
      elseif  solsta == MSK_SOL_STA_UNKNOWN
          println("Unknown solution status")
      else
          println("Other solution status")
      end

  end

  return primobj, dualobj
end
