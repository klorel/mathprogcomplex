using JuMP, SCS

m = Model(solver = SCSSolver())

vars_SDP = Dict{String, Array{JuMP.Variable, 2}}()
for i=1:4
    vars_SDP["toto$i"] = @variable(m, [1:3,1:3], SDP)
end

@objective(m, Max, 1 - sum( trace(transpose(A0) * Zi) for (cstrname, Zi) in vars_SDP))

@constraint(m, 1 - sum( trace(transpose(A1) * Zi) for (cstrname, Zi) in vars_SDP) == 0)

solve(m)



A0 = [2 -0.5 -0.6; -0.5 2 0.4; -0.6 0.4 3];
A1 = [0 1 0; 1 0 0; 0 0 0];
A2 = [0 0 1; 0 0 0; 1 0 0];
A3 = [0 0 0; 0 0 1; 0 1 0];

using JuMP, SCS
items  = [:Gold, :Silver, :Bronze]
values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)

m = Model(solver=SCSSolver())
@variable(m, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(m, Max, sum( values[item] * take[item] for item in items))
@constraint(m, sum( weight[item] * take[item] for item in items) <= 3)
solve(m)
println(getvalue(take))
# [Bronze] = 0.9999999496295456
# [  Gold] = 0.99999492720597
# [Silver] = 0.4666851698368782
