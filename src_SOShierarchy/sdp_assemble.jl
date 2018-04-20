function build_SDP(relaxctx::RelaxationContext, sdp::SDPInstance, mmtrelax_pb::MomentRelaxationPb)
    

    for ((cstrname, blocname), mmt) in mmtrelax_pb.constraints
        
        for ((γ, δ), poly) in mmt.mm
            for (expo, λ) in poly
                # Check the current monomial has correct degree
                if (relaxctx.hierarchykind==:Complex) && ((expo.degree.explvar > relaxctx.di[cstrname]) || (expo.degree.conjvar > relaxctx.di[cstrname]))
                    warn("convertMMtobase(): Found exponent pair of degree $(expo.degree) > $(relaxctx.di[cstrname]) for Complex hierarchy.\n($((α, β)), at $((γ, δ)) of MM matrix)")
                elseif (relaxctx.hierarchykind==:Real) && ((expo.degree.explvar > 2*relaxctx.di[cstrname]) || (expo.degree.conjvar != 0))
                    warn("convertMMtobase(): Found exponent pair of degree $(expo.degree) > 2*$(relaxctx.di[cstrname]) for Real hierarchy.\n($((α, β)), at $((γ, δ)) of MM matrix)")
                end
                !isnan(λ) || warn("convertMMtobase(): isNaN ! constraint $cstrname - clique $blocname - mm entry $((γ, δ)) - moment $((α, β))")

                # Determine which moment to affect the current coefficient.
                # NOTE: Hankel property could be exploited here, in real variables.
                α, β = split_expo(relaxctx, expo)

                # Add the current coeff to the SDP problem
                key = (cstrname, blocname, α, β)
                if !haskey(sdpbody, key)
                    sdpbody[key] = SortedDict{Tuple{Exponent, Exponent}, Number}()
                end
                Bi = sdpbody[key]

                if !haskey(Bi, (γ, δ))
                    Bi[(γ, δ)] = 0.0
                end
                Bi[(γ, δ)] += λ
            end
        end
    end

    for (expo, λ) in mmtrelax_pb.objective
        # Determine which moment to affect the current coefficient.
        # NOTE: Hankel property could be exploited here, in real variables.
        α, β = split_expo(relaxctx, expo)

        if !haskey(sdprhs, (α, β))
            sdprhs[(α, β)] = 0.0
        end

        sdprhs[(α, β)] += λ
    end

    return sdpbody, sdprhs
end
