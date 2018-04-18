function get_scen_vars(pt::Point, scenario::String)
    a = Iterators.filter(x->contains(x.name, scenario), keys(pt))
    return Point(SortedDict([k=>pt[k] for k in a]))
end

function get_volt_vars(pt::Point)
    a = Iterators.filter(x->contains(x.name, "VOLT"), keys(pt))
    return Point(SortedDict([k=>pt[k] for k in a]))
end

function get_bin_vars(pt::Point)
    a = Iterators.filter(x->contains(x.name, "Bin"), keys(pt))
    return Point(SortedDict([k=>pt[k] for k in a]))
end

function get_prod_vars(pt::Point)
    a = Iterators.filter(x->contains(x.name, "Sgen"), keys(pt))
    return Point(SortedDict([k=>pt[k] for k in a]))
end

function get_delta_var(pt::Point, scenario::String)
    a = Iterators.filter(x->contains(x.name, get_delta_varname(scenario)), keys(pt))
    return Point(SortedDict([k=>pt[k] for k in a]))
end

function get_splitted_Cpt(pt, scenario)
    vars = get_scen_vars(pt, scenario)
    volt_vars = real2cplx(get_volt_vars(vars))
    bin_vars = get_bin_vars(vars)
    prod_vars = real2cplx(get_prod_vars(vars))

    return volt_vars, bin_vars, prod_vars
end


###################
## Plot functions
###################

function plot_Volt_vars(point, scenarios)
    r1 = 0.9
    plot(r1*cos.(-π:0.01:π), r1*sin.(-π:0.01:π), label=:Vmin)
    r2 = 1.1
    plot!(r2*cos.(-π:0.01:π), r2*sin.(-π:0.01:π), label=:Vmax)
    for sc in scenarios
        println(sc)
        volt_vars, _, _ = get_splitted_Cpt(point, sc)
        ptsV = Array{Complex}(collect(values(volt_vars)))
        scatter!(real(ptsV), imag(ptsV), lab=sc)
    end
    gui()
end



function plot_Sgen_vars(pt_global, scenarios)
    fig = plot()
    pt = real2cplx(get_prod_vars(pt_global))
    Sgens = SortedDict()
    # Get prod by scenarios
    for scenario in scenarios
        _, _, bc_prod_vars = get_splitted_Cpt(pt_knitro, "BaseCase")
        Sgens[scenario] = bc_prod_vars
    end

    # One color per bus/generator
    init_scen = fisrt(Sgens).first
    i=1
    cols = SortedDict()
    for (varsc, valsc) in Sgens[init_scen]
        scenario, busid, genname, varname = split(varsc.name, "_")
        cols[(busid, genname)] = i
        i+=1
    end

    # Loop on generators
    for (varsc, valsc) in Sgens[init_scen]
        base_scenario, busid, genname, varname = split(varsc.name, "_")
        println("$base_scenario, $busid, $genname, $varname")

        gen = OPFpbs[base_scenario].ds.bus[bus_name(parse(busid))][genname]
        Smin, Smax = gen.power_min, gen.power_max
        domain = [ Smin, real(Smin)+im*imag(Smax), Smax, real(Smax)+im*imag(Smin), Smin]
        plot!(real(domain), imag(domain), color=cols[(busid, genname)], label="Cstr_$(busid)_$genname")

        for scenario in setdiff(scenarios, init_scen)
            gen_cur = OPFpbs[scenario].ds.bus[bus_name(parse(busid))][genname]
            Smin_cur, Smax_cur = gen_cur.power_min, gen_cur.power_max
            (Smin == Smin_cur) && (Smax == Smax_cur) || warn("bus $busid, gen $genname: between $base_scenario and $scenario, different bounds($Smin, $Smax vs. $Smin_cur, $Smax_cur)")

            var = Variable(variable_name("Sgen", bus_name(parse(busid)), String(genname), scenario), Complex)
            scatter!([pt[var]], color=cols[(busid, genname)], label=var.name)
        end
    end
    plot!(xlabel="P gen", ylabel = "Q gen")
    return fig
end


function plot_ViVj_vars(pt, scenarios)
    r1 = 0.9^2
    plt = plot(r1*cos.(-π:0.01:π), r1*sin.(-π:0.01:π), label=:Vmin2)
    r2 = 1.1^2
    plot!(r2*cos.(-π:0.01:π), r2*sin.(-π:0.01:π), label=:Vmax2)

    cols = SortedDict()
    i = 1
    for scen in scenarios
        cols[scen] = i
        i+=1
    end

    for scen in scenarios
        volt_vars, _, _ = get_splitted_Cpt(pt_knitro, scen)
        V = Vector{Complex}([val for (var, val) in volt_vars])
        ViVj = collect(SortedSet(V*V'))
        scatter!(ViVj, label=scen, color=cols[scen])
    end
    return plt
end
