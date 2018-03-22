struct SDPForm
    scal
    mats::Dict{String, AbstractMatrix}
end

function SDPForm(obj, B_i, expo)
    # println("----> SDPForm : $expo")
    scal = haskey(obj, expo) ? obj[expo] : 0
    mats = Dict{String, AbstractMatrix}()
    for (cstrname, mmb) in B_i
        if haskey(mmb.basis, expo)
            mats[cstrname] = mmb.basis[expo]
        end
    end

    # println("scal = $scal")
    # println("Matrices: $(keys(mats))")
    return SDPForm(scal, mats)
end


mutable struct SDPPrimal
    variables::Dict{String, Int}
    objective::SDPForm
    constraints::Dict{Exponent, SDPForm}
end

"""
    SDP_SOS = build_SDP_SOS(problem, di, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)

    Build the primal SDP corresponding to the dual SOS hierarchy for the provided problem.
"""
function build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)
    println("\n=== build_SDP_SOS(problem, max_cliques, B_i, cliquevarsbycstr, orderbyclique, relax_ctx)")
    println("Build the primal SDP corresponding to the dual SOS hierarchy for the provided problem.")

    # NOTE: vars should be read in args.
    d = relax_ctx.di["moment_cstr"]

    vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)


    cstrs = Dict{Exponent, SDPForm}()
    for realexpo in realexpos
        for conjexpo in conjexpos
            expo = product(realexpo, conjexpo)
            # println("==> Exponent $expo")

            cstr = SDPForm(problem.objective, B_i, expo)
            if (cstr.scal != 0) || length(cstr.mats) != 0
                cstrs[expo] = cstr
            else
                warn("Not adding $expo (no cstr...)")
            end
        end
    end

    obj = cstrs[Exponent()]
    delete!(cstrs, Exponent())

    vars = Dict{String, Int}()
    for (cstrname, mmb) in B_i
        vars[cstrname] = size(mmb.basis[Exponent()], 1)
    end
    SDPsizes = collect(values(vars))

    nbmatcstr = Int64[]
    for (expo, SDPform) in cstrs
        push!(nbmatcstr, length(SDPform.mats))
    end

    println("-> Nb of SDP variables:                                $(length(problem.constraints))")
    println("-> Size of the SDP variables:                          $(mean(SDPsizes)) / $(std(SDPsizes)) / $(maximum(SDPsizes)) (mean/std/max)")
    println("-> Nb of eq. SDP constraints:                          $(length(cstrs))")
    println("-> Nb of matrices by constraint:                       $(mean(nbmatcstr)) / $(std(nbmatcstr)) / $(maximum(nbmatcstr)) (mean/std/max)")
    println("-> Nb of coupling y_α,β variables:                     xx")
    println("-> Nb of constraint involved by coupling variables:    xx/xx (mean/std)")
    return SDPPrimal(vars, obj, cstrs)
end