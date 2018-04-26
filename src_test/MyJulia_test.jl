include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))

function read_args(ARGS)
    if length(ARGS)!=3
        error("First argument must be data_path for example ..\..\data\data_GOC
            Second argument must be instance folder for example Phase_0_IEEE14
            Third argument must be scenario for example scenario_1
        ")
    else
        data_path = ARGS[1]
        folder = ARGS[2]
        scenario = ARGS[3]
    end
    return data_path, folder, scenario
end

data_path, folder, scenario = read_args(ARGS)

function MyJulia(rawFile, genFile, contFile)
  println(rawFile, genFile, contFile)

  folder_path, scenario = splitdir(splitdir(rawFile)[1])
  folder = splitdir(folder_path)[2]

  outpath = joinpath("..","JuMP_runs","$(folder)_$(scenario)")
  ##read and load files
  OPFpbs = load_OPFproblems(rawFile, genFile, contFile)
  introduce_Sgenvariables!(OPFpbs)
  ## Building optimization problem
  pb_global = build_globalpb!(OPFpbs)

  ## convert into real problem
  pb_global_real = pb_cplx2real(pb_global)

  ##convert to JuMP model
  ##phase 1 : continuous relaxation
  println("
##############################################\n
Phase 1 : resolution  of continuous relaxation\n
##############################################")

  mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          # KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          # KTR_PARAM_FEASTOLABS=1.001e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=1)
  tic()
  my_timer = @elapsed m, variables_jump, ctr_jump, ctr_exp = get_JuMP_cartesian_model(pb_global_real, mysolver)
  @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  toc()
  #resolution
  solve(m)

  minslack = +Inf
  ctr_minslack = ""
  for (ctrname, (exp,lb,ub)) in ctr_exp
    body_value = getvalue(exp)
    slack = min(body_value-lb, ub-body_value)
    if slack <= minslack
      minslack = slack
      ctr_minslack = ctrname
    end
  end
  println("min slack for constraints in problem JuMP Phase 1 : ($minslack,$ctr_minslack)")

  ##phase 2 : resolution with complementary constraints + initial point = solution continuous relaxation
  println("
##############################################################################################\n
Phase 2 : resolution with complementary constraints from the solution of continuous relaxation\n
##############################################################################################")

  mysolver2 = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          # KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          # KTR_PARAM_FEASTOLABS=1.001e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=2)


  setsolver(m, mysolver2)
  # tic()
  # my_timer = @elapsed m2, variables_jump2 = get_JuMP_cartesian_model(pb_global_real, mysolver)
  # @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  # f = open(joinpath(outpath,"JuMP_solution_phase1.csv"),"w")
  # write(f, "Varname ; Value\n")
  # for (varname, varjump) in variables_jump
  #   write(f, "$varname; $(getvalue(varjump))\n")
  #   setvalue(variables_jump[varname], getvalue(varjump))
  # end
  # close(f)
  # toc()
  #resolution
  solve(m)

  ##get values
  println("Objective value : ", getobjectivevalue(m),"\n")

  # println("----Solution csv writing")
  # f = open(joinpath(outpath,"JuMP_solution.csv"),"w")
  # write(f, "Varname ; Value\n")
  # for (varname, var) in variables_jump
  #   value = getvalue(var)
  #   write(f, "$varname; $value\n")
  # end
  # close(f)

  minslack = +Inf
  ctr_minslack = ""
  # f = open(joinpath(outpath,"KnitroJuMP_constraints_eval.csv"),"w")
  # write(f, "Ctrname ; Value ; LB ; UB\n")
  for (ctrname, (exp,lb,ub)) in ctr_exp
    body_value = getvalue(exp)
    slack = min(body_value-lb, ub-body_value)
    if slack <= minslack
      minslack = slack
      ctr_minslack = ctrname
    end
    # write(f, "$ctrname; $body_value ; $lb ; $ub ; $slack\n")
  end
  # close(f)
  # println("----End solution csv writing\n")
  println("min slack for constraints in problem JuMP Phase 2 : ($minslack,$ctr_minslack)")

  ##create solution1.txt and solution2.txt
  println("--Solution txt writing")
  write_solutions(OPFpbs, variables_jump, outpath)
  println("--End solution txt writing\n")

  return pb_global_real

end



raw = "powersystem.raw"
gen = "generator.csv"
con = "contingency.csv"


instance_path = joinpath(pwd(),data_path, folder, scenario)



rawfile = joinpath(instance_path,raw)
genfile = joinpath(instance_path, gen)
contfile = joinpath(instance_path, con)


originalSTDOUT = STDOUT
date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
outpath = joinpath("..","JuMP_runs","$(folder)_$(scenario)")
isdir(outpath) || mkpath(outpath)
outlog = open(joinpath(outpath, "MyJulia1_$(date).log"), "w")
redirect_stdout(outlog)

println("Julia/JuMP test")

pb_global_real = MyJulia(rawfile,genfile,contfile)

# println("--Reading JuMP_solution.csv")
# sol = readdlm(joinpath(outpath,"JuMP_solution.csv"), ';',skipstart=1)
# pt = Point()
# for i in 1:size(sol,1)
#     varname = sol[i,1]
#     value = sol[i,2]
#     typevar = Real
#     if ismatch(r"Bin", varname)
#         typevar = Bool
#         if value > 0.5
#             value = 1.0
#         else
#             value = 1e-16
#         end
#     end
#     add_coord!(pt, Variable(varname,typevar) , value)
# end
#
#
# f = open(joinpath(outpath,"Polynomialpb_constraints_eval.csv"),"w")
# write(f, "Ctrname ; Value\n")
# for (ctrname, ctr) in pb_global_real.constraints
# body_value = evaluate(ctr.p, pt)
# write(f, "$ctrname; $body_value\n")
# end
# close(f)
#
# #
# # slacks = get_slacks(pb_global_real, pt)
# # for (ctr, s) in slacks
# #     if Real(s)<0
# #         println(ctr, " : ", s)
# #     end
# # end
# println("get_minslack solution point Julia polynomial problem :", get_minslack(pb_global_real, pt))

close(outlog)
redirect_stdout(originalSTDOUT)

println("--Reading solution1.txt and solution2.txt")
using MAT
sol_txt = read_solution_point_GOC(instance_path, outpath)
pt_txt = cplx2real(sol_txt)
# println(pt_txt[Variable("Scen1_6_Gen1_Sgen_Re",Real)])
# for (var, value) in pt_txt
#     if var.kind == Bool
#         println(var.name, value)
#     end
# end
# println(pt_txt == pt)

println("get_minslack point from txt files:", get_minslack(pb_global_real, pt_txt))
