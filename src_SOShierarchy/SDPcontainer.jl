function build_SDP(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxationPb)
    sdpbody = SDPBody()
    sdprhs = SDPRhs()

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


function print(io::IO, sdpbody::SDPBody)
    cstrlen = maximum(x->length(x[1]), keys(sdpbody))
    bloclen = maximum(x->length(x[2]), keys(sdpbody))
    expolen = maximum(x->length("($(x[3]), $(x[4]))"), keys(sdpbody))

    for ((cstrname, blocname, α, β), Bi) in sdpbody
        coordlen = maximum(x->length("($(x[1]), $(x[2]))"), keys(Bi))
        for ((γ, δ), λ) in Bi
            print_string(io, cstrname, cstrlen)
            print_string(io, blocname, bloclen)
            print_string(io, "($α, $β)", expolen); print(io, ": ")
            print_string(io, "($γ, $δ)", coordlen, alignright=false)
            println(io, "$λ")
        end
    end
end

function print(io::IO, sdprhs::SDPRhs)
    expolen = maximum(x->length("($(x[1]), $(x[2]))"), keys(sdprhs))
    for ((α, β), λ) in sdprhs
        print_string(io, "($α, $β)", expolen)
        println(io, "\t$λ")
    end
end

function get_maxlenkey(dict::SortedDict{String, U}) where U
    return maximum(map(x->length(x), keys(dict)))
end

function get_maxlenkey(dict::SortedDict{Tuple{Exponent,Exponent}, U}) where U
    return maximum(map(x->((α, β)=x; length("($α, $β)")), keys(dict)))
end