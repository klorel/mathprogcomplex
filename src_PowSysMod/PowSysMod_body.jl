include(joinpath("..", "src_PolynomialOptim", "PolynomialOptim.jl"))
using LightXML, MAT

import Base.copy

"""
    Link(orig, dest)

Represents a link from a bus origin `orig` to a bus destination `dest`
Util for stock or treat data about transmission lines or transformers
### Examples
```
julia > Link("BUS_1","BUS_4") #represents a link between BUS_1 and BUS_4
```
"""
struct Link
  orig::String
  dest::String
end

"""
    Abstract type AbstractNodeLabel

Abstract type for any element about a node
### Examples
```julia
struct MatpowerLoad <: AbstractNodeLabel
struct GOCVolt <: AbstractNodeLabel
struct IIDMGenerator <: AbstractNodeLabel
struct GOCShunt <: AbstractNodeLabel

```
"""
abstract type AbstractNodeLabel end

"""
    abstract type AbstractLinkLabel

Abstract type for any element about a link
### Examples
```
struct MatpowerLine_π <: AbstractLinkLabel
struct GOCLineπ_withtransformer <: AbstractLinkLabel
struct GOCNullImpedance_notransformer <: AbstractLinkLabel
```
"""
abstract type AbstractLinkLabel end

"""
    abstract type AbstractInput
Abstract type for any entry format
### Examples
```
type MatpowerInput <: AbstractInput end
type GOCInput <:AbstractInput end
type IIDMInput <: AbstractInput end
```
"""
abstract type AbstractInput end


"""
    DataSource(bus::Dict{String, Dict{String, Any}}, link::Dict{Link, Dict{String, Any}})

Store network data either in an attribute `bus` or in an attribute `link` depending on data type.

### Attributes
- `bus` : Dict{String, Dict{String, Any}}
- `link`: Dict{Link, Dict{String, Any}}

### Example
```
Instance Matpower WB2 with two nodes
julia > databus1 = Dict("Volt" => MatpowerVolt("BUS_1",1,0.95, 1.05), "Gen_1"=> MatpowerGenerator(-400 im, 600+400im, 3, [0 2 0], true))
julia > databus2 = Dict("Volt" => MatpowerVolt("BUS_2",2,0.95, 1.05), "Load" => MatpowerLoad(350-350im))
julia > busdata = Dict("BUS_1"=> databus1, "BUS_2" => databus2)
julia > link1to2 = Dict("LinkMP_1" => MatpowerLine_π(100,0.04,0.2,0.0,0.0,0.0))
julia > linkdata = Dict( Link("BUS_1","BUS_2") => link1to2)
julia > ds = DataSource(busdata,linkdata)
```
"""
struct DataSource
  # baseMVA::Number
  bus::Dict{String, Dict{String, Any}}
  link::Dict{Link, Dict{String, Any}}
end



"""
    struct GridStructure(scenario::String, node_linksin::Dict{String, Set{Link}}, node_linksout::Dict{String, Set{Link}})

Store information about network to construct power balances.

### Attributes
- `scenario::String` : scenario to study
- `generator_types::Set{Type}` : specify which types are generators
- `node_linksin::Dict{String, Set{Link}}` : for each node, set of links having this node for destination
- `node_linksout::Dict{String, Set{Link}}` : for each node, set of links having this node for origin

### Example
```
Instance Matpower WB2 with two nodes
julia > scenario = "BaseCase"
julia > node_linksin = Dict( "BUS_2" => Set{Link}([Link("BUS_1","BUS_2")]) )
julia > node_linksout = Dict( "BUS_1" => Set{Link}([Link("BUS_1, "BUS_2")]) )
julia > gs = GridStructure(scenario,generator_types, node_linksin, node_linksout)
```
"""
mutable struct GridStructure
  scenario::String
  node_linksin::Dict{String, Set{Link}}
  node_linksout::Dict{String, Set{Link}}
end


"""
    struct MathematicalProgramming(node_formulations::Dict{String, Dict{String, Symbol}}, link_formulations::Dict{String, Dict{String, Symbol}}, node_vars::Dict{String, Dict{String, Variable}}, link_vars::Dict{Link, Dict{String, Variable}})

Store information about mathematical formulation.
A node element or a link element can be formulated with a minimal number of variables or with variables representing all the quantities linked to the element.
This formulation is represented by symbols in the dictonaries `node_formulations` or `link_formulations` depending on the type of element.
Store created variables in `node_vars` or `link_vars`.

### Example
```
Instance Matpower WB2 with two nodes
julia > node_formulations =
julia > link_formulations = Dict( Link("BUS_1","BUS_2") => Dict("LinkMP_1" => :NbMinVar))
julia > node_vars =
julia > link_vars =

```
"""
mutable struct MathematicalProgramming
  node_formulations::Dict{String, Dict{String, Symbol}} # :NbMinVar, :NewVar, :None
  link_formulations::Dict{Link, Dict{String, Symbol}}
  node_vars::Dict{String, Dict{String, Variable}}
  link_vars::Dict{Link, Dict{String, Variable}}
end


"""
    struct Scenario(ds::DataSource, gs::GridStructure, mp::MathematicalProgramming)


### Example
```
Instance Matpower WB2 with two nodes
julia > basecase_scenario = Scenario(ds, gs, mp)
```
"""
mutable struct Scenario
  ds::DataSource
  gs::GridStructure
  mp::MathematicalProgramming
end

"""
    OPFproblems

### Examples
```

```
"""
const OPFProblems = Dict{String, Scenario}

# files to include for each entry type
include(joinpath("..","src_GOC","GOC_files.jl"))
include(joinpath("..","src_Matpower","Matpower_files.jl"))
include(joinpath("..","src_IIDM","IIDM_files.jl"))


# common files
include("add_generator_variables.jl")
include("build_globalpb.jl")
include("build_Problem.jl")
include("check_feasibility.jl")
include("default_behaviour.jl")
include("load_OPFproblems.jl")
include("naming_conventions.jl")
include("PowSysMod_accessors.jl")
include("PowSysMod_operators.jl")
include("utils_plots.jl")

 