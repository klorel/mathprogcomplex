# Default functions

## 1. Variables creation
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T
  return
end

function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T
  return
end

## 2. Power balance
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T
  return ["", Polynomial(), 0, 0]
end

function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T
  return Polynomial()
end

function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T
  return Polynomial()
end

## 3. Constraints creation
function constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T
  return SortedDict{String, Constraint}()
end

function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T
  return SortedDict{String, Constraint}()
end

## 4. Element cost
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T
  return 0
end

function cost(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{Link, Variable}, Snode::Polynomial, lb, ub) where T
  return 0
end
