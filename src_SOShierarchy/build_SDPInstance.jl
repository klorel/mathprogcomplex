function build_SDPInstance(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxationPb)
    sdpblocks = SDPBlocks()
    sdplin = SDPLin()
    sdpcst = SDPcst()
    block_to_vartype = SortedDict{String, Symbol}()

    ## Build blocks dict
    for ((cstrname, cliquename), mmt) in mmtrelax_pb.constraints
        block_name = get_blockname(cstrname, cliquename, mmtrelax_pb)
        block_to_vartype[block_name] = mmt.matrixkind

        for ((γ, δ), poly) in mmt.mm
            for (expo, λ) in poly
                # Check the current monomial has correct degree
                if (relaxctx.hierarchykind==:Complex) && ((expo.degree.explvar > relaxctx.di[cstrname]) || (expo.degree.conjvar > relaxctx.di[cstrname]))
                    warn("convertMMtobase(): Found exponent pair of degree $(expo.degree) > $(relaxctx.di[cstrname]) for Complex hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                elseif (relaxctx.hierarchykind==:Real) && ((expo.degree.explvar > 2*relaxctx.di[cstrname]) || (expo.degree.conjvar != 0))
                    warn("convertMMtobase(): Found exponent pair of degree $(expo.degree) > 2*$(relaxctx.di[cstrname]) for Real hierarchy.\n($(expo), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn("convertMMtobase(): isNaN ! constraint $cstrname - clique $blocname - mm entry $((γ, δ)) - moment $(expo)")

                # Determine which moment to affect the current coefficient.
                α, β = split_expo(relaxctx, expo)

                # Add the current coeff to the SDP problem
                key = ((α, β), block_name, γ, δ)

                @assert !haskey(sdpblocks, key)
                # Constraints are fα - ∑ Bi.Zi = 0
                sdpblocks[key] = -λ
            end
        end
    end

    # Collecting all monomial keys involving an overlapping variable
    overlap_vars = SortedSet(keys(mmtrelax_pb.vars_overlap))

    overlapping_expos_to_cliques = SortedDict{Exponent, SortedSet{String}}()

    for ((ctrmomentname, cliquename), mmtmat) in mmtrelax_pb.constraints
        if ctrmomentname == get_momentcstrname()
            warn("--> $ctrmomentname  $cliquename")
            for (key, poly) in mmtmat.mm
                println("$poly")
                for (expo, λ) in poly

                    expo_vars = SortedSet(keys(expo))
                    println("   $expo")
                    println("           $expo_vars")
                    println(intersect(expo_vars, overlap_vars))
                    if length(intersect(expo_vars, overlap_vars)) > 0
                        warn("Coucou")
                        info("Adding $expo , clique $cliquename")
                        haskey(overlapping_expos_to_cliques, cliquename) || (overlapping_expos_to_cliques[cliquename] = SortedSet{String}())
                        push!(overlapping_expos_to_cliques[expo], cliquename)
                    end
                end

            end

        end
    end

    @show overlapping_expos_to_cliques
    @show overlap_vars
    error()

    # Build linear dict
    for (var, cliques) in mmtrelax_pb.vars_overlap

        # Linking all clique instance of the var to the first clique var
        cliqueref = first(cliques)
        for clique in setdiff(cliques, SortedSet([cliqueref]))
            lagmult = get_ccmultvar(var, cliqueref, clique)
            sdplin[((Exponent(), Exponent(get_varinclique(var, cliqueref))), Exponent(lagmult))] = 1
            sdplin[((Exponent(), Exponent(get_varinclique(var, clique))), Exponent(lagmult))] = -1
        end

        # NOTE other strategies for linking variables could be implemented here.
    end

    ## Build constants dict
    for (expo, fαβ) in mmtrelax_pb.objective
        # Determine which moment to affect the current coefficient.
        α, β = split_expo(relaxctx, expo)

        if !haskey(sdpcst, (α, β))
            sdpcst[(α, β)] = 0.0
        end

        # Constraints are fα - ∑ Bi.Zi = 0
        sdpcst[(α, β)] += fαβ
    end

    return SDPInstance(block_to_vartype, sdpblocks, sdplin, sdpcst)
end


"""
    α, β = split_expo(expo::Exponent)

    Split the exponent into two exponents of conjugated and explicit variables in the complex case.
    Real case is not supported yet.
"""
function split_expo(relaxctx::RelaxationContext, expo::Exponent)
    α, β = Exponent(), Exponent()

    for (var, deg) in expo
        add_expod!(α, Exponent(SortedDict(var=>Degree(0, deg.conjvar))))
        add_expod!(β, Exponent(SortedDict(var=>Degree(deg.explvar, 0))))
    end

    if (relaxctx.hierarchykind == :Real) && (α.degree != Degree(0,0))
        error("split_expo(): Inconsistent degree $α, $β found for $(relaxctx.hierarchykind) hierarchy.")
    end
    return α, β
end


function print(io::IO, sdpinst::SDPInstance)
    println(io, " -- SDP Blocks:")
    print(io, sdpinst.blocks)
    println(io, " -- linear part:")
    length(sdpinst.lin) == 0 ? println("") : print(io, sdpinst.lin)
    println(io, " -- cc lin part:")
    for (k,v) in sdpinst.lin_cc
        println(io, "$k  $v")
    end
    println(io, " -- const part:")
    print(io, sdpinst.cst)
    println(io, " -- mat var types:")
    for (blockname, blocktype) in sdpinst.block_to_vartype
        println(io, "   $blockname  \t $blocktype")
    end
end

function print(io::IO, sdpblocks::SDPBlocks)
    cstrlen = maximum(x->length(format_string(x[1][1], x[1][2])), keys(sdpblocks))
    blocklen = maximum(x->length(x[2]), keys(sdpblocks))
    rowlen = maximum(x->length(format_string(x[3])), keys(sdpblocks))
    collen = maximum(x->length(format_string(x[4])), keys(sdpblocks))

    for (((α, β), blockname, γ, δ), λ) in sdpblocks
        print_string(io, format_string(α, β), cstrlen)
        print_string(io, blockname, blocklen)
        print_string(io, format_string(γ), rowlen)
        print_string(io, format_string(δ), collen)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end

function print(io::IO, sdplin::SDPLin)
    cstrlen = maximum(x->length(format_string(x[1][1], x[1][2])), keys(sdplin))
    varlen = maximum(x->length(format_string(x[2])), keys(sdplin))
    for (((α, β), var), λ) in sdplin
        print_string(io, format_string(α, β), cstrlen)
        print_string(io, format_string(var), varlen)
        @printf(io, "% .16e\n", λ)
    end
end

function print(io::IO, sdpcst::SDPcst)
    cstrlen = maximum(x->length(format_string(x[1], x[2])), keys(sdpcst))
    for ((α, β), λ) in sdpcst
        print_string(io, format_string(α, β), cstrlen)
        @printf(io, "% .16e\n", λ)
    end
end

function get_maxlenkey(dict::SortedDict{String, U}) where U
    return maximum(map(x->length(x), keys(dict)))
end

function get_maxlenkey(dict::SortedDict{Tuple{Exponent,Exponent}, U}) where U
    return maximum(map(x->((α, β)=x; length("($α, $β)")), keys(dict)))
end