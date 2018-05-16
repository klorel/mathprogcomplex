"""
    struct GOCLineπ_withtransformer <: AbstractLinkLabel

Structure descendant of AbstractLinkLabel

# Fields
- `link::Link`
- `id::String`
- `resistance::Float64`: coefficient rs
- `reactance::Float64`: coefficient xs
- `susceptance::Float64`: coefficient bc
- `transfo_ratio::Float64`: ratio τ
- `transfo_phase::Float64`: phase θ
- `power_magnitude_max::Float64`: Smax

"""
struct GOCLineπ_withtransformer <: AbstractLinkLabel
    orig_id::Int64
    dest_id::Int64
    id::Any
    resistance::Float64 #rs
    reactance::Float64 #xs
    susceptance::Float64 #bc
    transfo_ratio::Float64 #τ
    transfo_phase::Float64 #θ
    power_magnitude_max::Float64 #Smax
end

"""
    create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: GOCLineπ_withtransformer

Create voltage variables for origin and destination of `link` in `link_vars`.
Return nothing.
"""
function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: GOCLineπ_withtransformer
    origid, destid = link.orig, link.dest
    link_vars["Volt_orig"] = Variable(variable_name("VOLT", origid, "", scenario), Complex)
    link_vars["Volt_dest"] = Variable(variable_name("VOLT", destid, "", scenario), Complex)
    return
end


"""
    Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCLineπ_withtransformer

Return the polynomial power at the origin of the line depending quadratically on the voltages :
```julia
τ = element.transfo_ratio
θ = element.transfo_phase *pi/180
rs = element.resistance
xs = element.reactance
bc = element.susceptance
ys = 1/(rs+im*xs)
Yff = (ys+im*bc/2)/(τ^2)
Yft = -ys*1/(τ*exp(-im*θ))
Sorig = (Yff*link_variables["Volt_orig"] + Yft*link_variables["Volt_dest"]) * baseMVA
```
"""
function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCLineπ_withtransformer
    return get_baseMVA(link.orig) * link_vars["Volt_orig"] * conj(get_IorigGOC(element, link_vars))
end


"""
    Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCLineπ_withtransformer

Return the polynomial power at the destination of the line depending quadratically on the voltages :
```julia
τ = element.transfo_ratio
θ = element.transfo_phase *pi/180
rs = element.resistance
xs = element.reactance
bc = element.susceptance
ys = 1/(rs+im*xs)
Ytf = -ys*1/(τ*exp(im*θ))
Ytt = ys+im*bc/2
Sdest = (Ytf*link_variables["Volt_orig"] + Ytt*link_variables["Volt_dest"]) * baseMVA
```
"""
function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T<:GOCLineπ_withtransformer
   return get_baseMVA(link.dest) * link_vars["Volt_dest"] * conj(get_IdestGOC(element, link_vars))
end


## 3. Constraints creation
function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T<:GOCLineπ_withtransformer
    cstrs = SortedDict{String, Constraint}()

    if elem_formulation == :NbMinVar
        Sor = Sorig(element, link, elemid, elem_formulation, link_vars)
        Sde = Sdest(element, link, elemid, elem_formulation, link_vars)
        Smax = element.power_magnitude_max

        cstrs[get_Smax_orig_cstrname()] = (abs2(real(Sor)) + abs2(imag(Sor))) << Smax^2
        cstrs[get_Smax_dest_cstrname()] = (abs2(real(Sde)) + abs2(imag(Sde))) << Smax^2

        cstrs[get_Smax_orig_cstrname()].precond = :sqrt
        cstrs[get_Smax_dest_cstrname()].precond = :sqrt
    elseif elem_formulation == :NoCtr

    else
        error("Unknown formulation $elem_formulation for link $link, element $elemid.")
    end
    return cstrs
end

##Util functions
function get_IorigGOC(element::T, link_variables::SortedDict{String, Variable}) where T<: GOCLineπ_withtransformer
  τ = element.transfo_ratio
  θ = element.transfo_phase *pi/180
  rs = element.resistance
  xs = element.reactance
  bc = element.susceptance
  ys = 1/(rs+im*xs)
  Yff = (ys+im*bc/2)/(τ^2)
  Yft = -ys*1/(τ*exp(-im*θ))
  return Yff*link_variables["Volt_orig"] + Yft*link_variables["Volt_dest"]
end

function get_IdestGOC(element::T, link_variables::SortedDict{String, Variable}) where T<: GOCLineπ_withtransformer
  τ = element.transfo_ratio
  θ = element.transfo_phase *pi/180
  rs = element.resistance
  xs = element.reactance
  bc = element.susceptance
  ys = 1/(rs+im*xs)
  Ytf = -ys*1/(τ*exp(im*θ))
  Ytt = ys+im*bc/2
  return Ytf*link_variables["Volt_orig"] + Ytt*link_variables["Volt_dest"]
end
