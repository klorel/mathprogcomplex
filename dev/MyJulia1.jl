include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))


function MyJulia1(rawfile,genfile,contfile)
  ##read and load files
  OPFpbs = load_OPFproblems(rawfile, genfile, confile)
  introduce_Sgenvariables!(OPFpbs)
  ## Bulding optimization problem
  pb_global = build_globalpb!(OPFpbs)
  pb_global_real = pb_cplx2real(pb_global)

  ##convert to JuMP model


  ##create solution1.txt and solution2.txt
end
