"""
    Functions to copy structures (DataSource, GridStructure)
"""

function copy(ds::DataSource)
    bus = deepcopy(ds.bus)
    link = deepcopy(ds.link)
    return DataSource(bus, link)
end

function copy(gs::GridStructure)
    scenario = gs.scenario
    node_linksin = deepcopy(gs.node_linksin)
    node_linksout = deepcopy(gs.node_linksout)
    return GridStructure(scenario, node_linksin, node_linksout)
end
