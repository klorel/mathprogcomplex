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
    orig_id::Int64
    dest_id::Int64
    id::Any
    susceptance::Float64 #bc
    transfo_ratio::Float64 #τ
    transfo_phase::Float64 #θ
    power_magnitude_max::Float64 #Smax
end

"""
    create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: GOCNullImpedance_withtransformer

Create voltage variables and power variables for origin and destination of `link` in `link_vars`.
Return nothing.
"""
function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: GOCNullImpedance_withtransformer
    origid, destid = link.orig, link.dest
    link_vars["Volt_orig"] = Variable(variable_name("VOLT", origid, "", scenario), Complex)
    link_vars["Volt_dest"] = Variable(variable_name("VOLT", destid, "", scenario), Complex)

    link_vars[elemid*"_Sorig"] = Variable(variable_name("Sorig", origid, "", scenario), Complex)
    # link_vars[elemid*"_Sdest"] = Variable(variable_name("Sdest", destid, "", scenario), Complex)
    #elem_formulation modified ?
    return
end


"""
    Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCNullImpedance_withtransformer

Return power variable Sorig * baseMVA.
"""
function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCNullImpedance_withtransformer
    return get_baseMVA(link.orig)*link_vars[elemid*"_Sorig"]
end

"""
    Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCNullImpedance_withtransformer

Return power variable Sdest * baseMVA.
"""
function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCNullImpedance_withtransformer
   return get_baseMVA(link.dest)*(link_vars[elemid*"_Sorig"]- im * element.susceptance * abs2(link_vars["Volt_dest"]))
end


## 3. Constraints creation
function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <:GOCNullImpedance_withtransformer
    cstrs = SortedDict{String, Constraint}()

    #Smax constraints
    Sor = Sorig(element, link, elemid, elem_formulation, link_vars)
    Sde = Sdest(element, link, elemid, elem_formulation, link_vars)
    Smax = element.power_magnitude_max

    cstrs[get_Smax_orig_cstrname()] = (abs2(real(Sor)) + abs2(imag(Sor))) << Smax^2
    cstrs[get_Smax_dest_cstrname()] = (abs2(real(Sde)) + abs2(imag(Sde))) << Smax^2

    cstrs[get_Smax_orig_cstrname()].precond = :sqrt
    cstrs[get_Smax_dest_cstrname()].precond = :sqrt

    ## constraint linking Vorig and Vdest
    τ = element.transfo_ratio
    θ = element.transfo_phase
    cstrs[get_NullImpVolt_cstrname()] = Polynomial(link_vars["Volt_orig"] - link_vars["Volt_dest"] * 1/τ * exp(-im*θ))==0

    return cstrs
end
