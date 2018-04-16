function write_solutions(OPFpbs, variables_jump2, outpath)
    open(joinpath(outpath,"solution1.txt"),"w") do f
      write(f, "--generation dispatch\nbus id,unit id,pg(MW),qg(MVar) \n");
      for (busname, elems) in OPFpbs[basecase_scenario_name()].ds.bus
        for (elemname,element) in elems
                if typeof(element) == GOCGenerator
                  bus =  element.busid
                  gen = element.id
                  Pgen = getvalue(variables_jump2[variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Re"])
                  Qgen = getvalue(variables_jump2[variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Im"])
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
               V_re = getvalue(variables_jump2[variable_name("VOLT", busname, "", scenario)*"_Re"])
               V_im = getvalue(variables_jump2[variable_name("VOLT", busname, "", scenario)*"_Im"])
               V_mod = abs(V_re + V_im * im)
               V_theta = angle(V_re + V_im * im)*180/pi
               volt_values[(scenario_id,bus)] = (V_mod, V_theta)
             end
           end
         end
      else
        ##delta values
        scenario_id = String(matchall(r"\d+", scenario)[1])
        delta_values[scenario_id] = getvalue(variables_jump2[get_delta_varname(scenario)])
        ##volt values
         for (busname, elems) in OPFpb.ds.bus
           for (elemid, element) in elems
             if typeof(element) == GOCVolt
               bus = element.busid
               # bus = String(matchall(r"\d+", element.busname)[1])
               V_re = getvalue(variables_jump2[variable_name("VOLT", busname, "", scenario)*"_Re"])
               V_im = getvalue(variables_jump2[variable_name("VOLT", busname, "", scenario)*"_Im"])
               V_mod = abs(V_re + V_im * im)
               V_theta = angle(V_re + V_im * im)*180/pi
               volt_values[(scenario_id,bus)] = (V_mod, V_theta)
             elseif typeof(element) == GOCGenerator
               bus =  element.busid
               gen = element.id
               # Pgen = getvalue(variables_jump2[variable_name("Sgen", busname, elemid, scenario)*"_Re"])
               Qgen = getvalue(variables_jump2[variable_name("Sgen", busname, elemid, scenario)*"_Im"])
               # println("bus ",bus, " : ", element.power_min," << ", Pgen + im * Qgen," << ", element.power_max)
               Qgen_scen_values[(scenario_id,gen,bus, gen)] = Qgen
             end
           end
         end
      end

      for (link, elems) in OPFpb.ds.link
        orig = link.orig
        dest = link.dest
        Vo_re = getvalue(variables_jump2[variable_name("VOLT", orig, "", scenario)*"_Re"])
        Vo_im = getvalue(variables_jump2[variable_name("VOLT", orig, "", scenario)*"_Im"])
        Vd_re = getvalue(variables_jump2[variable_name("VOLT", dest, "", scenario)*"_Re"])
        Vd_im = getvalue(variables_jump2[variable_name("VOLT", dest, "", scenario)*"_Im"])

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

    open(joinpath(outpath,"solution2.txt"),"w") do f
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


  function write_solutions(OPFpbs::OPFProblems, pt::Point, outpath::String)
      open(joinpath(outpath,"solution1.txt"),"w") do f
        write(f, "--generation dispatch\nbus id,unit id,pg(MW),qg(MVar) \n");
        for (busname, elems) in OPFpbs[basecase_scenario_name()].ds.bus
          for (elemname,element) in elems
                  if typeof(element) == GOCGenerator
                    bus =  element.busid
                    gen = element.id
                    Pgen = pt[Variable(variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Re",Real)]
                    Qgen = pt[Variable(variable_name("Sgen", busname, elemname, basecase_scenario_name())*"_Im",Real)]
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
                 V_re = pt[Variable(variable_name("VOLT", busname, "", scenario)*"_Re",Real)]
                 V_im = pt[Variable(variable_name("VOLT", busname, "", scenario)*"_Im",Real)]
                 V_mod = abs(V_re + V_im * im)
                 V_theta = angle(V_re + V_im * im)*180/pi
                 volt_values[(scenario_id,bus)] = (V_mod, V_theta)
               end
             end
           end
        else
          ##delta values
          scenario_id = String(matchall(r"\d+", scenario)[1])
          delta_values[scenario_id] = pt[Variable(get_delta_varname(scenario),Real)]
          ##volt values
           for (busname, elems) in OPFpb.ds.bus
             for (elemid, element) in elems
               if typeof(element) == GOCVolt
                 bus = element.busid
                 # bus = String(matchall(r"\d+", element.busname)[1])
                 V_re = pt[Variable(variable_name("VOLT", busname, "", scenario)*"_Re",Real)]
                 V_im = pt[Variable(variable_name("VOLT", busname, "", scenario)*"_Im",Real)]
                 V_mod = abs(V_re + V_im * im)
                 V_theta = angle(V_re + V_im * im)*180/pi
                 volt_values[(scenario_id,bus)] = (V_mod, V_theta)
               elseif typeof(element) == GOCGenerator
                 bus =  element.busid
                 gen = element.id
                 # Pgen = getvalue(variables_jump2[variable_name("Sgen", busname, elemid, scenario)*"_Re"])
                 Qgen = pt[Variable(variable_name("Sgen", busname, elemid, scenario)*"_Im",Real)]
                 # println("bus ",bus, " : ", element.power_min," << ", Pgen + im * Qgen," << ", element.power_max)
                 Qgen_scen_values[(scenario_id,gen,bus, gen)] = Qgen
               end
             end
           end
        end

        for (link, elems) in OPFpb.ds.link
          orig = link.orig
          dest = link.dest
          Vo_re = pt[Variable(variable_name("VOLT", orig, "", scenario)*"_Re",Real)]
          Vo_im = pt[Variable(variable_name("VOLT", orig, "", scenario)*"_Im",Real)]
          Vd_re = pt[Variable(variable_name("VOLT", dest, "", scenario)*"_Re",Real)]
          Vd_im = pt[Variable(variable_name("VOLT", dest, "", scenario)*"_Im",Real)]

          link_elems_formulations = OPFpb.mp.link_formulations[link]
          link_vars = OPFpb.mp.link_vars[link]

          pt2 = Point()
          add_coord!(pt2, link_vars["Volt_orig"] , Vo_re + im * Vo_im)
          add_coord!(pt2, link_vars["Volt_dest"], Vd_re + im * Vd_im)

          for (elemid, element) in elems
              scenario == basecase_scenario_name() ? scenario_id = "0" : scenario_id = String(matchall(r"\d+", scenario)[1])
              link_id = element.id
              orig_id = element.orig_id
              dest_id = element.dest_id
              elem_formulation = link_elems_formulations[elemid]
              Sor = evaluate(Sorig(element, link, elemid, elem_formulation, link_vars),pt2)
              Sde = evaluate(Sdest(element, link, elemid, elem_formulation, link_vars),pt2)
              Slink_values[(scenario_id, link_id, orig_id, dest_id, link_id)] = (real(Sor),imag(Sor),real(Sde),imag(Sde))
           end
         end

       end

      open(joinpath(outpath,"solution2.txt"),"w") do f
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
