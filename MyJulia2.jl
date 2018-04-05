include(joinpath(pwd(),"src_PowSysMod", "PowSysMod_body.jl"))

function MyJulia2(rawFile, genFile, contFile)
  ##read and load files
  OPFpbs = load_OPFproblems(rawFile, genFile, contFile)
  introduce_Sgenvariables!(OPFpbs)
  ## Building optimization problem
  pb_global = build_globalpb!(OPFpbs)

  ## convert into real problem
  pb_global_real = pb_cplx2real(pb_global)

  ##convert to JuMP model
  mysolver = KnitroSolver(KTR_PARAM_OUTLEV=3,
                          KTR_PARAM_MAXIT=600,
                          KTR_PARAM_SCALE=0,
                          KTR_PARAM_FEASTOL=1.0,
                          KTR_PARAM_OPTTOL=1.0,
                          KTR_PARAM_FEASTOLABS=1e-6,
                          KTR_PARAM_OPTTOLABS=1e-3,
                          KTR_PARAM_BAR_INITPT=2,
                          KTR_PARAM_PRESOLVE=0,
                          KTR_PARAM_HONORBNDS=0,
                          KTR_PARAM_MIP_INTVAR_STRATEGY=2)
  tic()
  my_timer = @elapsed m, variables_jump = get_JuMP_cartesian_model(pb_global_real, mysolver)
  @printf("%-35s%10.6f s\n", "get_JuMP_cartesian_model", my_timer)
  toc()

  #resolution
  solve(m)

  ##get values
  println("Objective value : ", getobjectivevalue(m))

  # f = open("JuMP_solution.csv","w")
  # write(f, "Varname ; Value\n")
  # for (varname, var) in variables_jump
  #   value = getvalue(var)
  #   write(f, "$varname; $value\n")
  # end
  # close(f)

  ##create solution1.txt and solution2.txt
  println("Solution writing")

   open("solution1.txt","w") do f
     write(f, "--generation dispatch\nbus id,unit id,pg(MW),qg(MVar) \n");
     for (busname, elems) in OPFpbs[basecase_scenario_name()].ds.bus
       for (elemname,element) in elems
               if typeof(element) == GOCGenerator
                 bus =  element.busid
                 gen = element.id
                 # bus = matchall(r"\d+", element.busname)[1]
                 # gen = matchall(r"\d+", element.id)[1]
                 Pgen = getvalue(variables_jump[variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Re"])
                 Qgen = getvalue(variables_jump[variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Im"])
                 write(f, "$bus, $gen, $Pgen, $Qgen\n")
               end
        end
      end
     write(f,"--end of generation dispatch\n");
   end

   Qgen_scen_values = Dict{Tuple{String,Any,Int64,Any}, Float64}()
   volt_values = Dict{Tuple{String,Int64}, Tuple{Float64, Float64}}()
   delta_values = Dict{String,Float64}()
   Slink_values = Dict{Tuple{String,Any,Int64,Int64,Any}, Tuple{Float64, Float64, Float64, Float64}}()

   for (scenario, OPFpb) in OPFpbs
     if scenario==basecase_scenario_name()
       scenario_id = "0"
       ##volt values
        for (busname, elems) in OPFpb.ds.bus
          for (elemid, element) in elems
            if typeof(element) == GOCVolt
              bus = element.busid
              # bus = String(matchall(r"\d+", element.busname)[1])
              V_re = getvalue(variables_jump[variable_name("VOLT", busname, "", scenario)*"_Re"])
              V_im = getvalue(variables_jump[variable_name("VOLT", busname, "", scenario)*"_Im"])
              V_mod = abs(V_re + V_im * im)
              V_theta = angle(V_re + V_im * im)*180/pi
              volt_values[(scenario_id,bus)] = (V_mod, V_theta)
            end
          end
        end
     else
       ##delta values
       scenario_id = String(matchall(r"\d+", scenario)[1])
       delta_values[scenario_id] = getvalue(variables_jump[get_delta_varname(scenario)])
       ##volt values
        for (busname, elems) in OPFpb.ds.bus
          for (elemid, element) in elems
            if typeof(element) == GOCVolt
              bus = element.busid
              # bus = String(matchall(r"\d+", element.busname)[1])
              V_re = getvalue(variables_jump[variable_name("VOLT", busname, "", scenario)*"_Re"])
              V_im = getvalue(variables_jump[variable_name("VOLT", busname, "", scenario)*"_Im"])
              V_mod = abs(V_re + V_im * im)
              V_theta = angle(V_re + V_im * im)*180/pi
              volt_values[(scenario_id,bus)] = (V_mod, V_theta)
            elseif typeof(element) == GOCGenerator
              bus =  element.busid
              gen = element.id
               # bus = String(matchall(r"\d+", element.busname)[1])
               # gen = String(matchall(r"\d+", element.id)[1])
               Qgen = getvalue(variables_jump[variable_name("Sgen", busname, elemid, basecase_scenario_name())*"_Im"])
               Qgen_scen_values[(scenario_id,gen,bus, gen)] = Qgen
            end
          end
        end
     end

     for (link, elems) in OPFpb.ds.link
       orig = link.orig
       dest = link.dest
       Vo_re = getvalue(variables_jump[variable_name("VOLT", orig, "", scenario)*"_Re"])
       Vo_im = getvalue(variables_jump[variable_name("VOLT", orig, "", scenario)*"_Im"])
       Vd_re = getvalue(variables_jump[variable_name("VOLT", dest, "", scenario)*"_Re"])
       Vd_im = getvalue(variables_jump[variable_name("VOLT", dest, "", scenario)*"_Im"])

       link_elems_formulations = OPFpb.mp.link_formulations[link]
       link_vars = OPFpb.mp.link_vars[link]

       pt = Point()
       add_coord!(pt, link_vars["Volt_orig"] , Vo_re + im * Vo_im)
       add_coord!(pt, link_vars["Volt_dest"], Vd_re + im * Vd_im)

       for (elemid, element) in elems
           scenario == basecase_scenario_name() ? scenario_id = "0" : scenario_id = String(matchall(r"\d+", scenario)[1])
           link_id = element.id
           orig_id = element.orig_id
           dest_id = element.dest_id
           elem_formulation = link_elems_formulations[elemid]
           Sor = evaluate(Sorig(element, link, elemid, elem_formulation, link_vars),pt)
           Sde = evaluate(Sdest(element, link, elemid, elem_formulation, link_vars),pt)
           Slink_values[(scenario_id, link_id, orig_id, dest_id, link_id)] = (real(Sor),imag(Sor),real(Sde),imag(Sde))
        end
      end

    end

   open("solution2.txt","w") do f
     write(f, "--contingency generator\nconID,genID,busID,unitID,q(MW) \n");
     for ((sc, genID, busID, unitID),q) in Qgen_scen_values
       write(f,"$sc, $genID, $busID, $unitID, $q\n")
     end

     write(f,"--end of contingency generator\n--bus\ncontingency id,bus id,v(pu),theta(deg) \n");
     for ((sc,bus), (V_mod, V_theta)) in volt_values
       write(f, "$sc, $bus, $V_mod, $V_theta\n")
     end

     write(f,"--end of bus\n--Delta\ncontingency id,Delta(MW) \n");
     for (sc, delta) in delta_values
       write(f, "$sc, $delta \n")
     end

     write(f,"--end of Delta\n--line flow\ncontingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar) \n");
     for ((sc, br_id, orig_id, dest_id, br_id), (porig, qorig, pdest, qdest)) in Slink_values
       write(f, "$sc, $br_id, $orig_id, $dest_id, $br_id, $porig, $qorig, $pdest, $qdest \n")
     end
     write(f,"--end of line flow\n")
   end

end
