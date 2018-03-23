if length(ARGS)==0
  error("Argument must be a matpower instance as case9.m for example")
end

instance = ARGS[1]

instance_path = joinpath(pwd(),"..", "data_Matpower", "matpower",instance)

include(joinpath(pwd(),"..","src_PowSysMod", "PowSysMod_body.jl"))


function AMPLJulia_matpower(instance_path)
  ##read and load files
  OPFpbs = load_OPFproblems(MatpowerInput, instance_path)
  ## Bulding optimization problem
  pb_global = build_globalpb!(OPFpbs)

  tic()
  pb_global_real = pb_cplx2real(pb_global)
  t_cplx2real = toq()

  println("t_cplx2real : ", t_cplx2real)

  ## initial point
  init_point = Point()
  init_point_real = cplx2real(init_point)

  ## Exporting real problem
  tic()
  instance_name  = String(split(instance_path,'\\')[end])
  amplexportpath = joinpath("..","knitro_runs", "$(instance_name[1:end-2])")

  my_timer = @elapsed export_to_dat(pb_global_real, amplexportpath, init_point_real)
  @printf("%-35s%10.6f s\n", "export_to_dat", my_timer)
  t_buildexport = toq()

  println("time build export : ", t_buildexport)

  # _, t_knitro, _ = @timed run_knitro(amplexportpath, joinpath(pwd(),"..","src_ampl"))


end


AMPLJulia_matpower(instance_path)
