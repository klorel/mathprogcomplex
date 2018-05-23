function build_SDPInstance(relaxctx::RelaxationContext, mmtrelax_pb::MomentRelaxationPb)
    sdpblocks = SDPBlocks()
    sdplin = SDPLin()
    sdplinsym = SDPLinSym()
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
                # Constraints are fα - ∑ Bi.Zi = 0
                if mmt.matrixkind == :SDP || mmt.matrixkind == :CplxSDP
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
    cstrlenα= max(cstrlenα, length("#Ctr/Obj key j : conj part"))
    cstrlenβ = maximum(x->length(format_string(x[1][2])), keys(sdpblocks))
    cstrlenβ= max(cstrlenβ, length("#Ctr/Obj key j : expl part"))
    blocklen = maximum(x->length(x[2]), keys(sdpblocks))
    blocklen= max(blocklen, length("#Matrix variable key i"))
    rowlen = maximum(x->length(format_string(x[3])), keys(sdpblocks))
    rowlen = max(rowlen, length("#row key k"))
    collen = maximum(x->length(format_string(x[4])), keys(sdpblocks))
    collen = max(collen, length("#col key l"))

    print_string(io, "#Ctr/Obj key j : conj part", cstrlenα, indentedprint=indentedprint)
    print_string(io, "#Ctr/Obj key j : expl part", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "#Matrix variable key i", blocklen, indentedprint=indentedprint)
    print_string(io, "#row key k", rowlen, indentedprint=indentedprint)
    print_string(io, "#col key l", collen, indentedprint=indentedprint)
    @printf(io, "%23s %23s\n", "#A_ij[k, l] real part", "#A_ij[k, l] imag part")

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
    cstrlenα= max(cstrlenα, length("#Ctr/Obj key j : conj part"))
    cstrlenβ = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[1][2])), union(keys(sdplin), keys(sdplinsym))) : 0
    cstrlenβ= max(cstrlenβ, length("#Ctr/Obj key j : expl part"))

    varlen = length(sdplin)!=0 ? maximum(x->length(format_string(x[2])), keys(sdplin)) : 0
    varlensym = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[3], x[2])), keys(sdplinsym)) : 0
    varlen = max(varlen, varlensym, length("#Scalar variable key k"))

    print_string(io, "#Ctr/Obj key j : conj part", cstrlenα, indentedprint=indentedprint)
    print_string(io, "#Ctr/Obj key j : expl part", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "#Scalar variable key k", varlen, indentedprint=indentedprint)
    @printf(io, "%23s %23s\n", "#b_j[k] real part", "#b_j[k] imag part")

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


function print(io::IO, sdpcst::SDPcst; indentedprint=false)
    cstrlenα = maximum(x->length(format_string(x[1])), keys(sdpcst))
    cstrlenα= max(cstrlenα, length("#Ctr/Obj key j : conj part"))
    cstrlenβ = maximum(x->length(format_string(x[2])), keys(sdpcst))
    cstrlenβ= max(cstrlenβ, length("#Ctr/Obj key j : expl part"))

    print_string(io, "#Ctr/Obj key j : conj part", cstrlenα, indentedprint=indentedprint)
    print_string(io, "#Ctr/Obj key j : expl part", cstrlenβ, indentedprint=indentedprint)
    @printf(io, "%23s %23s\n", "#c_j real part", "#c_j imag part")

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