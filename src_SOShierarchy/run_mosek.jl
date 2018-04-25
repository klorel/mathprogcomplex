using Mosek

obj_key() = "1,1"

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
    if objctr == obj_key()
      nzc+=1
      barcj[nzc] = sdp_block.id
      barck[nzc] = upper
      barcl[nzc] = lower
      barcjkl[nzc] = coeff # * (lower==upper? 1: 0.5)
    else
      nza+=1
      barai[nza] = problem.name_to_ctr[objctr][1]
      baraj[nza] = sdp_block.id
      barak[nza] = upper
      baral[nza] = lower
      baraijkl[nza] = coeff # * (lower==upper? 1: 0.5)
    end
  end

  if debug
    warn("Reading SDP_Problem")
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


function solve_mosek(problem::SDP_Problem, primal::SortedDict{Tuple{String,String,String}, Float64}, dual::SortedDict{String, Float64}; debug = false)
  num_block = length(problem.id_to_block)
  # println("num_block : ", num_block)
  barvardim = [ length(problem.id_to_block[block].var_to_id) for block in 1:num_block ]
  vardim = 0
  # println(barvardim)
  println("num_block = ",   num_block)
  numcon, bkc, blc, buc = get_bounds(problem)
  println("numcon = ",   numcon)

  barcj, barck, barcl, barcjkl , barai, baraj, barak, baral, baraijkl = get_triplets(problem, debug=debug)
  # println(bkc)
  # println(blc)
  # println(buc)
  # Create a task object and attach log stream printer
  maketask() do task
      putstreamfunc(task,MSK_STREAM_LOG,printstream)

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
      putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

      # Set constraints matrices
      putbarablocktriplet(task, length(barai), barai, baraj, barak, baral, baraijkl)

      # Set contraints linear part
      # putaijlist(task, ai, aj, aij)

      # Objective matrices
      putbarcblocktriplet(task, length(barcj), barcj, barck, barcl, barcjkl)

      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE)
      # putintparam(task, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_AGGRESSIVE)

      if debug
        warn("Reading Mosek problem")
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
      end

      optimize(task)
      solutionsummary(task,MSK_STREAM_MSG)

      # analyzeproblem(task, MSK_STREAM_WRN)
      # analyzesolution(task, MSK_STREAM_LOG, MSK_SOL_ITR)
      # Get status information about the solution
      prosta = getprosta(task,MSK_SOL_ITR)
      solsta = getsolsta(task,MSK_SOL_ITR)
      if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
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

          println("Primal solution")
          for ((blockname, var1, var2), val) in primal
            @printf("%15s %5s %5s %f\n", blockname, var1, var2, val)
          end

          # Get dual solution
          println("Dual solution")
          for (ctrid, ctrname) in problem.id_to_ctr
            bars = getbarsj(task, MSK_SOL_ITR, ctrid)
            println("--> $ctrid, $ctrname")
            println(bars)
          end


          # activities = Dict{String, Float64}()
          # for ((objctr, block, var1, var2), coeff) in problem.matrices
          #   if !haskey(activities, objctr)
          #     activities[objctr] = 0
          #   end
          #   if haskey(primal, (block, var1, var2))
          #     activities[objctr] += coeff * primal[block, var1, var2]
          #   else
          #     activities[objctr] += coeff * primal[block, var2, var1]
          #   end
          #   # if objctr==obj_key()
          #   #   println(block, " - ", var1, " - ", var2, " - ", coeff)
          # end
          # # println(activities)
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
end
