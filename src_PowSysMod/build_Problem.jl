"""
    build_Problem!(OPFpbs, scenario::String)

Return a polynomial problem in complex variables with polynomial constraints
Step 0 : verify that there is at most one infered generator variable for each node
Step 1 : create variables in `OPFpbs[scenario].mp` according to formulations in `OPFpbs[scenario].mp` and add them in optimization problem `pb_opt`
Step 2 : create power balance constraints for each node in optimization problem `pb_opt` according to elements associated to the node
Step 3 : create objective in optimization problem `pb_opt` iterating on all nodal and link elements


# Arguments
- `OPFpbs` : dictionary of scenario (string) associated to structure Scenario
- `scenario` : String which must be a key of OPFpbs

# Output
- `pb_opt` : a polynomial optimization problem

# Example
```
julia > OPFpbs = build_OPFproblems(MatpowerInput, "instances\\matpower\\WB2.m")
julia > OPF = build_Problem!(OPFpbs, "BaseCase")
julia > print(OPF)
▶ variables: VOLT_1 VOLT_2
▶ objective:
(192.3) * conj(VOLT_1) * VOLT_1 + (-96.15384615384615 - 480.7692307692308im) * VOLT_1 * conj(VOLT_2) + (-96.15384615384615 + 480.7692307692308im) * conj(VOLT_1) * VOLT_2
▶ constraints:
LOAD_BUS_2: (-96.15384615384615 - 480.7692307692308im) * conj(VOLT_1) * VOLT_2 + (96.15384615384615 + 480.7692307692308im) * conj(VOLT_2) * VOLT_2 = -350.0f0 + 350.0f0im
UNIT_BUS_1: 0.0f0 - 400.0f0im < (96.15384615384615 + 480.7692307692308im) * conj(VOLT_1) * VOLT_1 + (-96.15384615384615 - 480.7692307692308im) * VOLT_1 * conj(VOLT_2) < 600.0f0 + 400.0f0im
   VOLTM_1: 0.9025 < conj(VOLT_1) * VOLT_1 < 1.1025
   VOLTM_2: 0.9025 < conj(VOLT_2) * VOLT_2 < 1.056784

```
"""


function build_Problem!(OPFpbs#=::OPFProblems=#, scenario::String)

  ds = OPFpbs[scenario].ds
  gs = OPFpbs[scenario].gs
  mp = OPFpbs[scenario].mp

  # 1. Creation des variables
  for bus in get_buses(OPFpbs, scenario)
    bus_elems = mp.node_formulations[bus]
    bus_vars = mp.node_vars[bus]
    for (elemid, elem) in ds.bus[bus]
      create_vars!(elem, bus, elemid, bus_elems[elemid], bus_vars, scenario)
    end
  end

  for link in get_links(OPFpbs, scenario)
    link_elems = mp.link_formulations[link]
    link_vars = mp.link_vars[link]
    for (elemid, elem) in ds.link[link]
      create_vars!(elem, link, elemid, link_elems[elemid], link_vars, scenario)
    end
  end

  # 2. Création des bilans
  pb_opt = Problem()
  for (node, vars) in mp.node_vars
    for (elem, var) in vars
      add_variable!(pb_opt, var)
    end
  end

  for (link, vars) in mp.link_vars
    for (elem, var) in vars
      add_variable!(pb_opt, var)
    end
  end

  for (busid, bus_elems) in ds.bus
    bus_elems_formulations = mp.node_formulations[busid]
    bus_elems_var = mp.node_vars[busid]

    cstrnames, Snode, lb, ub = get_Snodal(busid, bus_elems, bus_elems_formulations, bus_elems_var)
    S_balance = Snode + get_Slinks_in(busid, ds, gs, mp) + get_Slinks_out(busid, ds, gs, mp)

    add_constraint!(pb_opt, cstrname_nodal_balance(cstrnames, scenario, busid), lb << S_balance << ub)
  end


  ## 3. Création des contraintes
  for bus in get_buses(OPFpbs, scenario)
    bus_elems_formulations = mp.node_formulations[bus]
    bus_vars = mp.node_vars[bus]
    for (elemid, elem) in ds.bus[bus]
      elem_formulation = bus_elems_formulations[elemid]
      constraints = constraint(elem, bus, elemid, elem_formulation, bus_vars, scenario, OPFpbs)
      for (cstrname, cstr) in constraints
    	   add_constraint!(pb_opt, get_cstrname(scenario, bus, elemid, cstrname), cstr)
      end
    end
  end

  for link in get_links(OPFpbs, scenario)
    link_elems_formulations = mp.link_formulations[link]
    link_vars = mp.link_vars[link]
    for (elemid, elem) in ds.link[link]
      elem_formulation = link_elems_formulations[elemid]
      constraints = constraint(elem, link, elemid, elem_formulation, link_vars, scenario, OPFpbs)
      for (cstrname, cstr) in constraints
        add_constraint!(pb_opt, get_cstrname(scenario, link, elemid, cstrname), cstr)
      end
    end
  end

    ## 4. Création de l'objectif
    objective = Polynomial()
    for bus in get_buses(OPFpbs, scenario)
      bus_elems = ds.bus[bus]
      bus_elems_formulations = mp.node_formulations[bus]
      bus_elems_var = mp.node_vars[bus]

      _, Snode, lb, ub = get_Snodal(bus, bus_elems, bus_elems_formulations, bus_elems_var)
      add!(Snode, get_Slinks_in(bus, ds, gs, mp) + get_Slinks_out(bus, ds, gs, mp))

      for (elemid, elem) in bus_elems
        gencost::Polynomial = cost(elem, bus, elemid, bus_elems_formulations[elemid], bus_elems_var, Snode, lb, ub)
        add!(objective, gencost)
      end
    end

    set_objective!(pb_opt, objective)
  return pb_opt
end


"""
    get_Snodal(bus::String, bus_elems::Dict{String, Any}, bus_elems_formulations::Dict{String, Symbol}, bus_elems_var::Dict{String, Variable})

Return the total injected power at `bus` coming from nodal elements with bounds. The injected power for each specific element is computed using `Snodal`.

# Arguments
- `bus::String` :
- `bus_elems::Dict{String,Any}`:
- `bus_elems_formulations` :
- `bus_elems_var` :

# Output


#Examples
```
Instance Matpower WB2 with two nodes
julia > get_Snodal("BUS_1", )
```
"""
function get_Snodal(bus::String, bus_elems::Dict{String, Any}, bus_elems_formulations::Dict{String, Symbol}, bus_elems_var::Dict{String, Variable})
  cstrnames, Sres, lbres, ubres = [Set{String}(), Polynomial(), 0, 0]
  for (elemid, element) in bus_elems
    elem_formulation = bus_elems_formulations[elemid]
    cstrname_, Sres_, lbres_, ubres_ = Snodal(element, bus, elemid, elem_formulation, bus_elems_var)
    Sres, lbres, ubres = [Sres, lbres, ubres] + [Sres_, lbres_, ubres_]
    push!(cstrnames, cstrname_)
  end

  # Remove all constraints names equal to ""
  dense_cstrnames = Set{String}()
  for cstr in cstrnames
    if cstr != ""
      push!(dense_cstrnames, cstr)
    end
  end
  return [dense_cstrnames, Sres, lbres, ubres]
end

function get_Slinks_in(bus, ds::DataSource, gs::GridStructure, mp::MathematicalProgramming)
  Sres = Polynomial()
  for link in get_links_in(bus, gs)
    link_elems_formulation = mp.link_formulations[link]
    link_vars = mp.link_vars[link]
    for (elemid, elem) in ds.link[link]
      elem_formulation = link_elems_formulation[elemid]
      add!(Sres, Sdest(elem, link, elemid, elem_formulation, link_vars))
    end
  end
  return Sres
end

function get_Slinks_out(bus, ds::DataSource, gs::GridStructure, mp::MathematicalProgramming)
  Sres = Polynomial()
  for link in get_links_out(bus, gs)
    link_elems_formulations = mp.link_formulations[link]
    link_vars = mp.link_vars[link]
    for (elemid, elem) in ds.link[link]
      elem_formulation = link_elems_formulations[elemid]
      add!(Sres, Sorig(elem, link, elemid, elem_formulation, link_vars))
    end
  end
  return Sres
end
