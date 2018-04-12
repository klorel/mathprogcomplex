## Add a dimension for the cliques in moment_matrices...
function build_SDP(relaxctx::RelaxationContext, moment_matrices::SortedDict{Tuple{String, String}, MomentMatrix})
    sdpbody = SDPBody()
    sdprhs = SDPRhs()

    for ((cstrname, cliquename), mm) in moment_matrices
        
        for ((γ, δ), poly) in mm
            for (expo, λ) in poly
                # Determine which moment to affect the current coefficient.
                # NOTE: could symmetry properties be exploited here ? Is this the same thing ?
                α, β = split_expo(relaxctx, expo)

                # Check the current monomial has correct degree
                if (β.degree.explvar > relctx.di[cstrname]) || (α.degree.conjvar > relctx.di[cstrname])
                    warn("convertMMtobase(): Found exponent pair of degree $(α.degree), $(β.degree) > $(relctx.di[cstrname]) ($((α, β)), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn("convertMMtobase(): isNaN ! constraint $cstrname - clique $cliquename - mm entry $((γ, δ)) - moment $((α, β))")

                add_to_dict!(sdpbody[(cstrname, cliquename, α, β)], (γ, δ), λ)
            end
        end
    end

    for (expo, λ) in pb.objective
        # Determine which moment to affect the current coefficient.
        # NOTE: could symmetry properties be exploited here ? Is this the same thing ?
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
            add_expod!(α, var^expo.conjdeg)
            add_expod!(β, var^expo.expldeg)
        end
    else
        ## TODO: C'est ici que ça se joue en réel...
        error("build_SDP(): Real hierarchy not yet handled.")
    end
    
    return α, β
end