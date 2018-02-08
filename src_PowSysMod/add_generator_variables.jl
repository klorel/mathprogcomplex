"""
    introduce_Sgenvariables!(OPFpbs)

Modify `mp` field of all `Scenario` in `OPFpbs` in order to create Sgen variables when the problem will be constructed.
"""
function introduce_Sgenvariables!(OPFpbs)
    for (scenario, OPFpb) in OPFpbs
        node_formulations = OPFpbs[scenario].mp.node_formulations

        for (bus, elems_form) in node_formulations
            for elem in keys(elems_form)
                if ismatch(r"Gen", elem)
                    if scenario != basecase_scenario_name()
                        elems_form[elem] = :GOCcoupling
                    else
                        elems_form[elem] = :NewVar
                    end
                elseif ismatch(r"Load", elem) || ismatch(r"Shunt", elem)
                    elems_form[elem] = :NbMinVar
                else
                    elems_form[elem] = :NewVar
                end
            end
        end
    end
    return OPFpbs
end
