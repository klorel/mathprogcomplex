function build_SDPInstance(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxationPb)
    sdpblocks = SDPBlocks()
    sdplin = SDPLin()
    sdplinsym = SDPLinSym()
    sdpcst = SDPCst()
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
                # Constraints are fα - ∑ Bi.Zi = 0
                if mmt.matrixkind == :SDP || mmt.matrixkind == :SDPC
                    key = ((α, β), block_name, γ, δ)
                    @assert !haskey(sdpblocks, key)

                    sdpblocks[key] = -λ
                elseif mmt.matrixkind == :Sym || mmt.matrixkind == :CplxSym
                    key = ((α, β), block_name, product(γ, δ))
                    haskey(sdplinsym, key) || (sdplinsym[key] = 0)

                    sdplinsym[key] += -λ * (γ!=δ ? 2 : 1)
                else
                    error("build_SDPInstance(): Unhandled matrix kind $(mmt.matrixkind) for ($cstrname, $cliquename)")
                end

            end
        end
    end

    ## Build linear dict
    ## TODO : enforce clique coupling constraints

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

    return SDPInstance(block_to_vartype, sdpblocks, sdplinsym, sdplin, sdpcst)
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
    println(io, " -- const part:")
    print(io, sdpinst.cst)
    println(io, " -- mat var types:")
    for (blockname, blocktype) in sdpinst.block_to_vartype
        println(io, "   $blockname  \t $blocktype")
    end
end

function print(io::IO, sdpblocks::SDPBlocks; indentedprint=false)
    cstrlenα = maximum(x->length(format_string(x[1][1])), keys(sdpblocks))
    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ = maximum(x->length(format_string(x[1][2])), keys(sdpblocks))
    cstrlenβ= max(cstrlenβ, length("j_expl"))
    blocklen = maximum(x->length(x[2]), keys(sdpblocks))
    blocklen= max(blocklen, length("Zi"))
    rowlen = maximum(x->length(format_string(x[3])), keys(sdpblocks))
    rowlen = max(rowlen, length("k"))
    collen = maximum(x->length(format_string(x[4])), keys(sdpblocks))
    collen = max(collen, length("l"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "Zi", blocklen, indentedprint=indentedprint)
    print_string(io, "k", rowlen, indentedprint=indentedprint)
    print_string(io, "l", collen, indentedprint=indentedprint)
    print_string(io, "Real(A_ij[k, l])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(A_ij[k, l])", 23, indentedprint=indentedprint)
    println(io)

    for (((α, β), blockname, γ, δ), λ) in sdpblocks
        print_string(io, format_string(α), cstrlenα, indentedprint=indentedprint)
        print_string(io, format_string(β), cstrlenβ, indentedprint=indentedprint)
        print_string(io, blockname, blocklen, indentedprint=indentedprint)
        print_string(io, format_string(γ), rowlen, indentedprint=indentedprint)
        print_string(io, format_string(δ), collen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end

function print(io::IO, sdplinsym::SDPLinSym, sdplin::SDPLin; indentedprint=false)
    cstrlenα = length(sdplin)!=0 ? maximum(x->length(format_string(x[1][1])), union(keys(sdplin), keys(sdplinsym))) : 0
    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[1][2])), union(keys(sdplin), keys(sdplinsym))) : 0
    cstrlenβ= max(cstrlenβ, length("j_expl"))

    varlen = length(sdplin)!=0 ? maximum(x->length(format_string(x[2])), keys(sdplin)) : 0
    varlensym = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[3], x[2])), keys(sdplinsym)) : 0
    varlen = max(varlen, varlensym, length("x[k]"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "x[k]", varlen, indentedprint=indentedprint)
    print_string(io, "Real(b_j[k])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(b_j[k])", 23, indentedprint=indentedprint)
    println(io)

    if length(sdplin)!=0
        for (((α, β), var), λ) in sdplin
            print_string(io, format_string(α), cstrlenα, indentedprint=indentedprint)
            print_string(io, format_string(β), cstrlenβ, indentedprint=indentedprint)
            print_string(io, format_string(var), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
        end
    end
    if length(sdplinsym) != 0
        for (((α, β), blockname, var), λ) in sdplinsym
            print_string(io, format_string(α), cstrlenα, indentedprint=indentedprint)
            print_string(io, format_string(β), cstrlenβ, indentedprint=indentedprint)
            print_string(io, format_string(var, blockname), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
        end
    end
end


function print(io::IO, sdpcst::SDPCst; indentedprint=false)
    cstrlenα = maximum(x->length(format_string(x[1])), keys(sdpcst))
    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ = maximum(x->length(format_string(x[2])), keys(sdpcst))
    cstrlenβ= max(cstrlenβ, length("j_expl"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "Real(c_j)", 23, indentedprint=indentedprint)
    print_string(io, "Imag(c_j)", 23, indentedprint=indentedprint)
    println(io)

    for ((α, β), λ) in sdpcst
        print_string(io, format_string(α), cstrlenα, indentedprint=indentedprint)
        print_string(io, format_string(β), cstrlenβ, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end

function get_maxlenkey(dict::SortedDict{String, U}) where U
    return maximum(map(x->length(x), keys(dict)))
end

function get_maxlenkey(dict::SortedDict{Tuple{Exponent,Exponent}, U}) where U
    return maximum(map(x->((α, β)=x; length("($α, $β)")), keys(dict)))
end