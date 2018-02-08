"""
    Functions to copy structures (DataSource, GridStructure)
"""

function copy(ds::DataSource)
    bus = copy(ds.bus)
    link = copy(ds.link)
    return DataSource(bus, link)
end

function copy(gs::GridStructure)
    scenario = gs.scenario
    node_linksin = copy(gs.node_linksin)
    node_linksout = copy(gs.node_linksout)
    return GridStructure(scenario, node_linksin, node_linksout)
end
