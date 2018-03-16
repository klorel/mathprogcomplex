"""
    B_i_dict = compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    Compute the decomposition of localizing matrix corresponding to each constraint on the moment variable basis, yielding several matrices B_i,α,β for each constraint i.
"""
function compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)
    println("\n=== compute_Bibycstr(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)")
    println("Compute the decomposition of localizing matrix corresponding to each constraint on the moment variable basis, yielding several matrices B_i,α,β for each constraint i.")

    di = relax_ctx.di
    ki = relax_ctx.ki

    Bi_tot = Dict{String, MomentMatrixBasis}()
    # NOTE: should chose variables in the cliquevarsbycstr obj...
    vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])

    for (cstrname, cstr) in problem.constraints
        d, k = di[cstrname], ki[cstrname]
        mm = MomentMatrix(vars, d - k)

        Bi_tot[cstrname] = convertMMtobase(cstr.p * mm, d, k)
    end


    nb_matrices = Int64[]
    sparsity = Float64[]
    for (cstrname, mmb) in Bi_tot
        push!(nb_matrices, length(mmb.basis))
        for (expo, mat) in mmb.basis
            push!(sparsity, count(mat .!= 0))
        end
    end

    println("-> Nb matrices:                        $(sum(nb_matrices))")
    println("-> Nb matrices by constraint:          $(mean(nb_matrices)) / $(std(nb_matrices)) (mean/std)")
    println("-> Nb nnz by matrix:                   $(mean(sparsity)) / $(std(sparsity)) (mean/std)")
    return Bi_tot
end
