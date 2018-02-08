"""
    struct GOCNullImpedance_withtransformer <: AbstractLinkLabel

Structure descendant of AbstractLinkLabel

# Fields
- `link::Link`
- `id::String`
- `susceptance::Float64`: coefficient bc
- `transfo_ratio::Float64`: ratio τ
- `transfo_phase::Float64`: phase θ
- `power_magnitude_max::Float64`: Smax

"""
struct GOCNullImpedance_withtransformer <: AbstractLinkLabel
    link::Link
    id::String
    susceptance::Float64 #bc
    transfo_ratio::Float64 #τ
    transfo_phase::Float64 #θ
    power_magnitude_max::Float64 #Smax
end

"""
    create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}, scenario::String) where T <: GOCNullImpedance_withtransformer

Create voltage variables and power variables for origin and destination of `link` in `link_vars`.
Return nothing.
"""
function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}, scenario::String) where T <: GOCNullImpedance_withtransformer
    origid, destid = link.orig, link.dest
    link_vars["Volt_orig"] = Variable(variable_name("VOLT", origid, "", scenario), Complex)
    link_vars["Volt_dest"] = Variable(variable_name("VOLT", destid, "", scenario), Complex)

    link_vars[elemid*"_Sorig"] = Variable(variable_name("Sorig", origid, "", scenario), Complex)
    link_vars[elemid*"_Sdest"] = Variable(variable_name("Sdest", destid, "", scenario), Complex)
    #elem_formulation modified ?
    return
end


"""
    Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:GOCNullImpedance_withtransformer

Return power variable Sorig * baseMVA.
"""
function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:GOCNullImpedance_withtransformer
    return get_baseMVA(link.orig)*link_vars[elemid*"_Sorig"]
end

"""
    Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:GOCNullImpedance_withtransformer

Return power variable Sdest * baseMVA.
"""
function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}) where T<:GOCNullImpedance_withtransformer
   return get_baseMVA(link.dest)*link_vars[elemid*"_Sdest"]
end


## 3. Constraints creation
# function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::Dict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MyType
#   return Dict{String, Constraint}()
# end
