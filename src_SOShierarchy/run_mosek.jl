using Mosek

obj_key() = "1,1"

printstream(msg::String) = print(msg)

function get_triplets(problem::SDP_Problem)

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
      baraijkl[nza] = coeff * (lower==upper? 1: 0.5)
    end
  end
  # println( "barcj : ", barcj)
  # println( "barck : ", barck)
  # println( "barcl : ", barcl)
  # println( "barcjkl : ", barcjkl)
  #
  # println( "barai : ", barai)
  # println( "baraj : ", baraj)
  # println( "barak : ", barak)
  # println( "baral : ", baral)
  # println( "baraijkl : ", baraijkl)

  return barcj, barck, barcl, barcjkl, barai, baraj, barak, baral, baraijkl
end
function get_bounds(problem::SDP_Problem)
  numcon=length(problem.name_to_ctr)
  MOSEK_KIND = Dict(["EQ"=>MSK_BK_FX, "GEQ"=>MSK_BK_LO, "LEQ"=>MSK_BK_UP, "RNG"=>MSK_BK_RA])
  bkc = Boundkey[ Mosek.Boundkey(1)  for kv in problem.name_to_ctr]
  # constraint
  buc = Float64[0 for i in 1:numcon]
  blc = Float64[0 for i in 1:numcon]
  for kv in problem.name_to_ctr
    # @printf("%10d%20s%10d%20s\n", kv[2][1], kv[1], MOSEK_KIND[kv[2][2]], kv[2][2])
    id_ctr=kv[2][1]
    lb = kv[2][3]
    ub = kv[2][4]
    cst = 0
    if haskey(problem.cst_ctr, kv[1])
      cst = problem.cst_ctr[kv[1]]
    end    
    bkc[id_ctr] = MOSEK_KIND[kv[2][2]]
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
      error("get_bounds() : Unknown constraint kind $(kv[2][2]) $(bkc[id_ctr]) $(MSK_BK_FX[1])")
    end
    # if bkc[id_ctr]!= MSK_BK_UP
    #   blc[id_ctr] = rhs
    # end
    # if bkc[id_ctr]!= MSK_BK_LO
    #   buc[id_ctr] = rhs
    # end
  end
  return numcon, bkc, blc, buc
end


function solve_mosek(problem::SDP_Problem, primal::Dict{Tuple{String,String,String}, Float64}, dual::Dict{String, Float64})
  num_block = length(problem.id_to_block)
  # println("num_block : ", num_block)
  barvardim = [ length(problem.id_to_block[block].var_to_id) for block in 1:num_block ]
  vardim = 0
  # println(barvardim)
  println("num_block = ",   num_block)
  numcon, bkc, blc, buc = get_bounds(problem)
  println("numcon = ",   numcon)

  barcj, barck, barcl, barcjkl , barai, baraj, barak, baral, baraijkl = get_triplets(problem)
  # println(bkc)
  # println(blc)
  # println(buc)
  # Create a task object and attach log stream printer
  maketask() do task
      putstreamfunc(task,MSK_STREAM_LOG,printstream)
  
      println("EE")

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

      optimize(task)
      solutionsummary(task,MSK_STREAM_MSG)

      # analyzeproblem(task, MSK_STREAM_WRN)
      # analyzesolution(task, MSK_STREAM_LOG, MSK_SOL_ITR)
      # Get status information about the solution
      prosta = getprosta(task,MSK_SOL_ITR)
      solsta = getsolsta(task,MSK_SOL_ITR)
      if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
          # Output a solution
          for id_block in problem.id_to_block
            barx = getbarxj(task, MSK_SOL_ITR, id_block[1])
            all_variables = ["" for kv in id_block[2].var_to_id]
            for kv in id_block[2].var_to_id
              all_variables[kv[2]]=kv[1]
            end
            n = 0
            for j in 1:length(all_variables)
              for i in j:length(all_variables)
                n+=1
                # @printf("%15s%15s%15s%25.10f\n", id_block[2].name, all_variables[i], all_variables[j], barx[n])
                primal[id_block[2].name, all_variables[i], all_variables[j]] = barx[n]
              end
            end
            # primal[block.name, ]
          end
          # @printf("Optimal solution: \n  xx = %s\n  barx = %s\n", xx',barx')
          activities = Dict{String, Float64}()
          for ((objctr, block, var1, var2), coeff) in problem.matrices
          # for objctr_block_var1_var2_coeff in problem.matrices
            # objctr = objctr_block_var1_var2_coeff[1]
            if !haskey(activities, objctr)
              activities[objctr] = 0
            end
            # for block_var1_var2_coeff in objctr_block_var1_var2_coeff[2]
              # block = block_var1_var2_coeff[1][1]
              # var1 = block_var1_var2_coeff[1][2]
              # var2 = block_var1_var2_coeff[1][3]
              # coeff = block_var1_var2_coeff[2]
              if haskey(primal, (block, var1, var2))
                activities[objctr] += coeff * primal[block, var1, var2]
              else
                activities[objctr] += coeff * primal[block, var2, var1]
              end
              # if objctr==obj_key()
              #   println(block, " - ", var1, " - ", var2, " - ", coeff)
              # end
            # end
          end
          # println(activities)
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
