include(joinpath(ROOT,"..","src_PowSysMod", "PowSysMod_body.jl"))

if length(ARGS) != 2
  error("Expected input:\n  julia compare_dats.jl file1.dat file2.dat")
end

file1, file2 = ARGS

maxerror, errors = compare_dat(file1, file2, 1e-10)

count = 0
for key in keys(errors)
  count += length(errors[key])
end

println("Max error on coeff: $maxerror")
println("Nb of coeffs with error > 1e-10: $(count)")
