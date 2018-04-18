# Expanding ComplexModeler

### Reading the new files(s)
In `/src/ComplexModeler/read_functions/`, a function `read_input(intput_type::T, filenames::Array{String})` allowing to build the datastructure holding all relevant numerical values has to be implemented. It must return :
- baseMVA::Float64: a float value allowing to convert the power quantities from p.u. (per unit) to MVA (MegaVolt Ampere),
- bus::SortedDict{String, MyType1}(): a dictionnary associating a bus name to the corresponding numerical values (load, generator, ...), stored in a custom datatype MyType1,
- link::SortedDict{Link, MyType2}(): a dictionnary associating a link to the corresponding numerical values (resistance, reactance, ...), stored in a custom data type,
- node_elems::SortedDict{String, SortedSet{Any}}(): a dictionnary associating each bus to a set of "labels" that describes it. Each "label" is a custom julia type, corresponding to a nodal element such as generator, load, shunt, ...
- link_elems::SortedDict{Link, SortedSet{Any}}(): a dictionnary associating each link to a set of "labels" that describes it. Each "label" is a custom julia type, corresponding to a nodal element such as Line_π, NullImpedance, ...

A template including the function prototype, return arguments, and location of custom data storage type is provided.

### Accessing the data
In order to ease access to the stored data, accessors can be implemented in `/src/ComplexModeler/base_accessors/`. It is likely these functions will receive as arguments the bus or link, the OPFProblem structure holding all the data and the current scenario, e.g.
```julia
function get_Sload(bus::String, OPFpbs::OPFProblems, scenario::String)
```

### Building variables, power balance and constraints by element
To each nodal or link element introduced by the new data format, such as null impedance lines, new generator types, etc corresponds a new julia type, derived from `AbstractNodeLabel` or `AbstractLinkLabel` --- these types are used in node_elems and link_elems to describe the network. To each of these types must be associated a number of functions allowing to determine:
- their contribution to each nodal power balance (if any)
- what new variables they require
- the supplemental constraints (if any)

These functions have to be implemented in `/src/ComplexModeler/structural_fcts/`, a template is also provided (`template_elements.jl`).

**TODO: provide explicit API for adding data in GS (new variables)**

```julia
type MyType <: AbstractNodeLabel end

function create_vars!(elem::T, bus::String, OPFpbs::OPFProblems, scenario::String) where T <: Type{MyType}
  return
end
```
