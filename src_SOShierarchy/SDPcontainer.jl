function build_SDP(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxationPb)
    sdpbody = SDPBody()
    sdprhs = SDPRhs()

    for ((cstrname, cliquename), mmt) in mmtrelax_pb.constraints
        
        for ((γ, δ), poly) in mmt.mm
            for (expo, λ) in poly
                # Determine which moment to affect the current coefficient.
                # NOTE: Hankel property could be exploited here, in real variables.
                println("$(typeof(expo)) $expo  -> $λ")
                α, β = split_expo(relaxctx, expo)

                # Check the current monomial has correct degree
                if (β.degree.explvar > relaxctx.di[cstrname]) || (α.degree.conjvar > relaxctx.di[cstrname])
                    warn("convertMMtobase(): Found exponent pair of degree $(α.degree), $(β.degree) > $(relctx.di[cstrname]) ($((α, β)), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn("convertMMtobase(): isNaN ! constraint $cstrname - clique $cliquename - mm entry $((γ, δ)) - moment $((α, β))")

                key1 = (cstrname, cliquename, α, β)
                key2 = (γ, δ)
                if !haskey(sdpbody, key1)
                    sdpbody[key1] = Dict{Tuple{Exponent, Exponent}, Number}()
                end
                add_to_dict!(sdpbody[key1], (γ, δ), λ)
            end
        end
    end

    for (expo, λ) in mmtrelax_pb.objective
        # Determine which moment to affect the current coefficient.
        # NOTE: Hankel property could be exploited here, in real variables.
        α, β = split_expo(relaxctx, expo)
        add_to_dict!(sdprhs, (α, β), λ)
    end

    return sdpbody, sdprhs
end


"""
    α, β = split_expo(expo::Exponent)

    Split the exponent into two exponents of conjugated and explicit variables in the complex case.
    Real case is not supported yet.
"""
function split_expo(relaxctx::RelaxationContext, expo::Exponent)
    α, β = Exponent(), Exponent()

    if relaxctx.hierarchykind == :Complex
        for (var, deg) in expo
            add_expod!(α, Exponent(Dict(var=>Degree(0, expo.degree.conjvar))))
            add_expod!(β, Exponent(Dict(var=>Degree(expo.degree.explvar, 0))))
        end
    else
        ## TODO: C'est ici que ça se joue en réel...
        error("build_SDP(): Real hierarchy not yet handled.")
    end
    
    return α, β
end