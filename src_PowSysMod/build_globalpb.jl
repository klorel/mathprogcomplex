"""
     build_globalpb!(OPFpbs)

Return a polynomial problem in complex variables with polynomial constraints combining all scenarios contained in `OPFpbs`

# Arguments
- `OPFpbs` : dictionary of scenario (string)

# Output
- `pb_opt` : a polynomial optimization problem

# Example
```
julia > OPFpbs = build_OPFproblems(MatpowerInput, "instances\\matpower\\WB2.m")
julia > pb = build_globalpb!(OPFpbs)
```
"""

function build_globalpb!(OPFpbs)
      constraints = SortedDict{String,SortedDict{String,Constraint}}()
      variables = SortedDict{String,SortedDict{String, Type}}()
      pb_global = Problem()
      pb_global = build_Problem!(OPFpbs, basecase_scenario_name())
      for scenario in setdiff(collect(keys(OPFpbs)), [basecase_scenario_name()])
            pb_scenario = build_Problem!(OPFpbs, scenario)
            constraints[scenario] = pb_scenario.constraints
            variables[scenario] = pb_scenario.variables
      end
      for (scenario, vars) in variables
            for var in vars
                  add_variable!(pb_global, var)
            end
      end
      for (scenario, constraints) in constraints
            for (ctr_name, ctr) in constraints
                  add_constraint!(pb_global, ctr_name, ctr)
            end
      end
      return pb_global
end
