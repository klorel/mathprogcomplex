"""
    File containing all the conventions used to name objects.
"""
bus_name(id_bus::Int64) = "BUS_$id_bus"
get_busid(busname::String) = busname[5:end]
basecase_scenario_name() = "BaseCase"

function variable_name(varname::String, bus::String, elemid::String, scenario::String)
    if elemid == ""
        return "$(scenario)_$(get_busid(bus))_$(varname)"
    else
        return "$(scenario)_$(get_busid(bus))_$(elemid)_$(varname)"
    end
end

function variable_name(varname::String, link::Link, elemid::String, scenario::String)
    if elemid == ""
        return "$(scenario)_$(get_busid(link.orig))-$(get_busid(link.dest))_$(varname)"
    else
        return "$(scenario)_$(get_busid(link.orig))-$(get_busid(link.dest))_$(elemid)_$(varname)"
    end
end

function cstrname_nodal_balance(cstrnames::SortedSet{String}, scenario::String, bus::String)
    cstr = ""
    for cstrname in sort(collect(cstrnames))
        cstr = "$(cstr)-$cstrname"
    end

    if cstr == ""
      cstr = "-LOAD"
    end

    return "$(scenario)_$(get_busid(bus))_BALANCE$(cstr)"
end

function get_cstrname(scenario::String, bus::String, elemid::String, cstrname::String)
    if elemid == ""
        return "$(scenario)_$(get_busid(bus))_$(cstrname)"
    else
        return "$(scenario)_$(get_busid(bus))_$(elemid)_$(cstrname)"
    end
end

function get_cstrname(scenario::String, link::Link, elemid::String, cstrname::String)
    if elemid == ""
        return "$(scenario)_$(get_busid(link.orig))-$(get_busid(link.dest))_$(cstrname)"
    else
        return "$(scenario)_$(get_busid(link.orig))-$(get_busid(link.dest))_$(elemid)_$(cstrname)"
    end
end


### GOC conventions
scenarioname(id_contingency::Int) = "Scen$id_contingency"

# nodal element names
volt_name() = "Volt"
generator_name(gen_id) = "Gen$gen_id"
load_name(id_load) = "Load$id_load"
shunt_name(id_shunt) = "Shunt$id_shunt"

get_delta_varname(scenario::String) = "$(scenario)_Delta"
get_binInf_varname(scenario1::String, scenario2::String, bus::String) = "BinVolt_$(get_busid(bus))_$(scenario1)_inf_$(scenario2)"
get_binEq_varname(scenario1::String, scenario2::String, bus::String) = "BinVolt_$(get_busid(bus))_$(scenario1)_eq_$(scenario2)"

# constraint names
get_VoltM_cstrname() = "VOLTM"
get_ShuntDef_cstrname() = "ShuntDef"
get_LoadDef_cstrname() = "LoadDef"
get_GenBounds_cstrname() = "GenBounds"
get_NullImpVolt_cstrname() = "NullImpVolt"


get_VoltBinDef_upper() = "BinDef-upper"
get_VoltBinDef_lower() = "BinDef-lower"
get_VoltBinDef_complement() = "BinDef-compl"

get_CC_active_cstrname() = "CC-Pgen"
get_CC_reactiveupper_cstrname() = "CC-Qgen-upper"
get_CC_reactivelower_cstrname() = "CC-Qgen-lower"

get_Smax_orig_cstrname() = "Smax-orig"
get_Smax_dest_cstrname() = "Smax-dest"


get_GOC_Volt_ϵ() = 1e-6
get_GOC_BigM() = 1e1

##node elems

##link elems
nullimpedance_notransformer_name(branch_id) = "NullImp_notransfo_$branch_id"
linepi_withtransformer_name(branch_id) = "Lineπ_transfo_$branch_id"
nullimpedance_withtransformer_name(branch_id) = "NullImp_transfo_$branch_id"
linepi_notransformer_name(branch_id) = "Lineπ_notransfo_$branch_id"
