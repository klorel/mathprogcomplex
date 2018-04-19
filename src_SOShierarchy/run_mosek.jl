using Mosek
printstream(msg::String) = print(msg)
function get_triplets(problem::SDP_Problem)

  nzc = 0
  nza = 0
  for objctr_block_var1_var2_coeff in problem.matrix
    objctr = objctr_block_var1_var2_coeff[1]
    for block_var1_var2_coeff in objctr_block_var1_var2_coeff[2]
      var1 = block_var1_var2_coeff[1][2]
      var2 = block_var1_var2_coeff[1][3]

      @assert(var1!="NONE" && var2!="NONE")
      if objctr == "OBJ"
        nzc += 1
      else
        nza += 1
      end
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
  for objctr_block_var1_var2_coeff in problem.matrix
    objctr = objctr_block_var1_var2_coeff[1]
    for block_var1_var2_coeff in objctr_block_var1_var2_coeff[2]

      block = block_var1_var2_coeff[1][1]
      var1 = block_var1_var2_coeff[1][2]
      var2 = block_var1_var2_coeff[1][3]
      coeff = block_var1_var2_coeff[2]

      sdp_block = problem.name_to_block[block]
      lower = min(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
      upper = max(sdp_block.var_to_id[var1], sdp_block.var_to_id[var2])
      if objctr == "OBJ"
        nzc+=1
        barcj[nzc] = sdp_block.id
        barck[nzc] = upper
        barcl[nzc] = lower
        barcjkl[nzc] = coeff * (lower==upper? 1: 0.5)
      else
        nza+=1
        barai[nza] = problem.name_to_ctr[objctr][1]
        baraj[nza] = sdp_block.id
        barak[nza] = upper
        baral[nza] = lower
        baraijkl[nza] = coeff * (lower==upper? 1: 0.5)
      end
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
  numcon=length(sdp.name_to_ctr)
  MOSEK_KIND = Dict(["EQ"=>MSK_BK_FX, "GEQ"=>MSK_BK_LO, "LEQ"=>MSK_BK_UP, "RNG"=>MSK_BK_RA])
  bkc = Int32[ -1  for kv in sdp.name_to_ctr]
  # constraint
  buc = Float64[0 for i in 1:numcon]
  blc = Float64[0 for i in 1:numcon]
  for kv in sdp.name_to_ctr
    # @printf("%10d%20s%10d%20s\n", kv[2][1], kv[1], MOSEK_KIND[kv[2][2]], kv[2][2])
    id_ctr=kv[2][1]
    lb = kv[2][3]
    ub = kv[2][4]
    cst = sdp.cst_ctr[kv[1]]
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
      @printf("unknow constraint kind\n")
      exit(0)
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
      # Append matrix variables of sizes in 'BARVARDIM'.
      # The variables will initially be fixed at zero.
      appendbarvars(task,barvardim)
      # Append 'numcon' empty constraints.
      # The constraints will initially have no bounds.
      appendcons(task,numcon)
      # Set the bounds on constraints.
      putconboundslice(task,1,numcon+1, bkc,blc,buc)
      # Minimize
      putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
      #
      putbarablocktriplet(task, length(barai), barai, baraj, barak, baral, baraijkl)

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
          for objctr_block_var1_var2_coeff in problem.matrix
            objctr = objctr_block_var1_var2_coeff[1]
            if !haskey(activities, objctr)
              activities[objctr] = 0
            end
            for block_var1_var2_coeff in objctr_block_var1_var2_coeff[2]
              block = block_var1_var2_coeff[1][1]
              var1 = block_var1_var2_coeff[1][2]
              var2 = block_var1_var2_coeff[1][3]
              coeff = block_var1_var2_coeff[2]
              if haskey(primal, (block, var1, var2))
                activities[objctr] += coeff * primal[block, var1, var2]
              else
                activities[objctr] += coeff * primal[block, var2, var1]
              end
              # if objctr=="OBJ"
              #   println(block, " - ", var1, " - ", var2, " - ", coeff)
              # end
            end
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
#
# # Bound keys for constraints
# bkc = [MSK_BK_FX,
#        MSK_BK_FX]
#
# # Bound values for constraints
# blc = [1.0, 0.5]
# buc = [1.0, 0.5]
#
#
# A = sparse( [1,2,2],[1,2,3],[1.0, 1.0, 1.0])
# conesub = [1, 2, 3]
#
# barci = [1, 2, 2, 3, 3]
# barcj = [1, 1, 2, 2, 3]
# barcval = [2.0, 1.0, 2.0, 1.0, 2.0]
#
# barai   = Any[ [1, 2, 3],
#                [1, 2, 3, 2, 3, 3] ]
# baraj   = Any[ [1, 2, 3],
#                [1, 1, 1, 2, 2, 3] ]
# baraval = Any[ [1.0, 1.0, 1.0],
#                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] ]
#
# numvar = 3
# numcon = length(bkc)
# barvardim = [3]
#
# # Create a task object and attach log stream printer
# maketask() do task
#     putstreamfunc(task,MSK_STREAM_LOG,printstream)
#
#     # Append 'numvar' variables.
#     # The variables will initially be fixed at zero (x=0).
#     appendvars(task,numvar)
#
#     # Append 'numcon' empty constraints.
#     # The constraints will initially have no bounds.
#     appendcons(task,numcon)
#
#     # Append matrix variables of sizes in 'BARVARDIM'.
#     # The variables will initially be fixed at zero.
#     appendbarvars(task,barvardim)
#
#     # Set the linear term c_0 in the objective.
#     putcj(task, 1, 1.0)
#
#     # Set the bounds on variable j
#     # blx[j] <= x_j <= bux[j]
#     putvarboundslice(task,1,numvar+1,
#                      [ MSK_BK_FR::Int32 for i in 1:numvar ],
#                      [ -Inf             for i in 1:numvar ],
#                      [ +Inf             for i in 1:numvar ])
#
#     # Set the bounds on constraints.
#     # blc[i] <= constraint_i <= buc[i]
#     putconboundslice(task,1,numcon+1, bkc,blc,buc)
#
#     # Input row i of A
#     putacolslice(task,1,numvar+1,
#                  A.colptr[1:numvar], A.colptr[2:numvar+1],
#                  A.rowval,A.nzval)
#
#     appendcone(task,MSK_CT_QUAD, 0.0, conesub)
#
#     symc  = appendsparsesymmat(task,barvardim[1],
#                                barci,
#                                barcj,
#                                barcval)
#
#     syma0 = appendsparsesymmat(task,barvardim[1],
#                                barai[1],
#                                baraj[1],
#                                baraval[1])
#
#     syma1 = appendsparsesymmat(task,barvardim[1],
#                                barai[2],
#                                baraj[2],
#                                baraval[2])
#
#     putbarcj(task,1, [symc], [1.0])
#
#     putbaraij(task,1, 1, [syma0], [1.0])
#     putbaraij(task,2, 1, [syma1], [1.0])
#
#     # Input the objective sense (minimize/maximize)
#     putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
#
#     # Solve the problem and print summary
#     optimize(task)
#     solutionsummary(task,MSK_STREAM_MSG)
#
#     # Get status information about the solution
#     prosta = getprosta(task,MSK_SOL_ITR)
#     solsta = getsolsta(task,MSK_SOL_ITR)
#
#
#     if solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL
#         # Output a solution
#         xx = getxx(task,MSK_SOL_ITR)
#         barx = getbarxj(task,MSK_SOL_ITR, 1)
#
#         @printf("Optimal solution: \n  xx = %s\n  barx = %s\n", xx',barx')
#     elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
#         println("Primal or dual infeasibility.\n")
#     elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
#         println("Primal or dual infeasibility.\n")
#     elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
#         println("Primal or dual infeasibility.\n")
#     elseif  solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
#         println("Primal or dual infeasibility.\n")
#     elseif  solsta == MSK_SOL_STA_UNKNOWN
#         println("Unknown solution status")
#     else
#         println("Other solution status")
#     end
# end
