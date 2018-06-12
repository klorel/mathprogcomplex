###############################################################################
####  Relaxation context
###############################################################################
function print_build_relctx(relax_ctx, pb)
    outstream = []
    relaxparams = relax_ctx.relaxparams
    relaxparams[:opt_outmode]!=1 && push!(outstream, STDOUT)
    relaxparams[:opt_outmode]≥0  && push!(outstream, open(relaxparams[:opt_outname], "a"))

    (relaxparams[:opt_outlev] == 0) && return

    # Compute indicators
    nb_densecstrs = 0
    maxdeg_densecstr = Float32[]
    for (cstr, ki_) in relax_ctx.ki
        if relax_ctx.di[cstr] > ki_
            nb_densecstrs += 1
            push!(maxdeg_densecstr, ki_)
        end
    end

    exposet = Set()
    nb_expotot = 0
    degbycstr = Int64[]
    for (cstrname, cstr) in get_constraints(pb)
        push!(degbycstr, max(cstr.p.degree.explvar, cstr.p.degree.conjvar))
        for (expo, λ) in cstr.p
            push!(exposet, expo)
            nb_expotot += 1
        end
    end

    # print
    for outstr in outstream
        if relaxparams[:opt_outlev] ≥ 1
            println(outstr, "\n=== set_relaxation()")

            print(outstr, "Relaxation is ", string(relaxparams[:opt_hierarchykind]))
            relaxparams[:opt_issparse] && print(outstr, ", sparse")
            relaxparams[:opt_multiordered] && print(outstr, ", multiordered")
            println(outstr, ".")

            println(outstr, "Global order is : ", relaxparams[:opt_globalorder], ".")

            print(outstr, "PhaseInvariance ")
            !relaxparams[:opt_sym_phaseinv] && print(outstr, "not ")
            print(outstr, "required and ")
            !relaxparams[:pb_isphaseinv] && print(outstr, "not ")
            println(outstr, "found.")

            println(outstr, "\n=> Problem with:")
            println(outstr, "-> Nb of variables:                        ", length(pb.variables))
            println(outstr, "-> Nb of constraints:                      ", length(pb.constraints))
            println(outstr, "-> Nb of monomials:                        $nb_expotot ($(length(exposet)) different)")
            @printf(outstr, "-> Max total degree by constraint:         %.1f / %.1f (mean/std)\n", mean(degbycstr), std(degbycstr))

            println(outstr, "=> Relaxation characteristics:")
            @printf(outstr, "-> Number of constraints s.t. di > ki:     %i / %i\n", nb_densecstrs, length(pb.constraints))
            @printf(outstr, "-> Max total degree on such constraints:   %.1f / %.1f (mean/std)\n", mean(maxdeg_densecstr), std(maxdeg_densecstr))
            @printf(outstr, "All variables appearing in such constraints will be linked in the sparsity pattern, which will largely densify it.\n")
        end

        (outstr!=STDOUT) && close(outstr)
    end
end


function print(io::IO, relctx::RelaxationContext)
    print(io, "RelaxationContext:\n")
    print(io, "ismultiordered         : $(relctx.ismultiordered)\n")
    print(io, "issparse               : $(relctx.issparse)\n")
    print(io, "symmetries             : $(relctx.symmetries)\n")
    print(io, "hierarchykind          : $(relctx.hierarchykind)\n")
    print(io, "renamevars             : $(relctx.renamevars)\n")
    for cstrname in sort(collect(keys(relctx.di)))
        di = relctx.di[cstrname]
        print(io, "di                     : $cstrname  \t=> $di\n")
    end
    for cstrname in sort(collect(keys(relctx.ki)))
        ki = relctx.ki[cstrname]
        print(io, "ki                     : $cstrname  \t=> $ki\n")
    end
    for cstrname in sort(collect(keys(relctx.cstrtypes)))
        cstrtype = relctx.cstrtypes[cstrname]
        print(io, "bar var types          : $cstrname  \t=> $(string(cstrtype))\n")
    end
end


###############################################################################
####  Moment problem construction
###############################################################################
function print_build_momentrelax(relax_ctx, momentrelaxation, nb_expos)
    (relax_ctx[:opt_outlev] == 0) && return

    outstream = IOStream[]
    relax_ctx[:opt_outmode]!=1 && push!(outstream, STDOUT)
    relax_ctx[:opt_outmode]≥0  && push!(outstream, open(relax_ctx[:opt_outname], "a"))

    ## Compute indicators
    nb_overlap_expos = length(expo_to_cliques)

    for outstr in outstream
        if relax_ctx[:opt_outlev] ≥ 1
            println(outstr, "\n=== MomentRelaxation(relax_ctx, problem, moment_param::Dict{String, Tuple{Set{String}, Int}}, max_cliques::Dict{String, Set{Variable}})")
            println(outstr, "Compute the moment and localizing matrices associated with the problem constraints and clique decomposition and return a MomentRelaxation object.")

            if relax_ctx[:opt_outlev] ≥ 2
                print(outstr, momentrelax)
            end
            if nb_overlap_expos > 0
                info(outstr, "Nb exponents coupled: $nb_overlap_expos (over $nb_expos)")
            end

            ## NOTE: which relevant indicators here ?
        end

        (outstr!=STDOUT) && close(outstr)
    end
end


function print(io::IO, momentrelax::MomentRelaxation{T}) where T
    println(io, "Moment Relaxation Problem:")
    println(io, "→ Objective: ")
    momentlen = maximum(x->length(string(x)), keys(momentrelax.objective))
    for moment in sort(collect(keys(momentrelax.objective)))
        coeff = momentrelax.objective[moment]
        print_string(io, string(moment), momentlen)
        println(io, " $coeff")
    end

    println(io, "→ Constraints:")
    for (cstrname, blocname) in sort(collect(keys(momentrelax.constraints)))
        mmtmat = momentrelax.constraints[(cstrname, blocname)]
        println(io, " → $cstrname, $blocname")
        println(io, mmtmat)
    end

    println(io, "→ Moments clique overlap:")
    if length(momentrelax.moments_overlap) > 0
        mmtlength = maximum(x->length(string(x)), keys(momentrelax.moments_overlap))
        for moment in sort(collect(keys(momentrelax.moments_overlap)))
            cliquenames = momentrelax.moments_overlap[moment]
            print(io, " → ")
            print_string(io, string(moment), mmtlength)
            for clique in sort(collect(cliquenames)) print(io, "$clique, ") end
            @printf(io, "\b\b \n")
        end
    else
        print(io, "  None")
    end
end