function export_SDP(relax_ctx, sdp::SDPInstance, path; indentedprint=true, renamemoments=true)

    # Build moment shortname dict if required
    momentdict = build_momentdict(sdp, renamemoments)

    # Collect all constraints keys
    ctr_keys = build_ctrkeysset(sdp)

    # Export blocks of constraints
    blocks_file = joinpath(path, "matrix.sdp")
    !isfile(blocks_file) || rm(blocks_file)

    fblocks = open(blocks_file, "a")
    print_blocksfile(fblocks, sdp.blocks, momentdict, indentedprint=indentedprint)
    close(fblocks)

    # Export linear part of constraints
    lin_file = joinpath(path, "lin.sdp")
    !isfile(lin_file) || rm(lin_file)

    flin = open(lin_file, "a")
    print_linfile(flin, sdp.lin, sdp.linsym, momentdict, indentedprint=indentedprint)
    close(flin)

    # Export constants of constraints
    cst_file = joinpath(path, "const.sdp")
    !isfile(cst_file) || rm(cst_file)

    fcst = open(cst_file, "a")
    print_cstfile(fcst, sdp.cst, momentdict, ctr_keys, indentedprint=indentedprint)
    close(fcst)


    # Export bloc types
    types_file = joinpath(path, "types.sdp")
    !isfile(types_file) || rm(types_file)

    ftypes = open(types_file, "a")
    print_typesfile(ftypes, sdp.block_to_vartype)
    close(ftypes)

    # Export moment dictionary name if required
    momname_file = joinpath(path, "names.sdp")
    !isfile(momname_file) || rm(momname_file)

    if renamemoments
        fname = open(momname_file, "a")
        print_namesfile(fname, momentdict)
        close(fname)
    end
end


"""
    momentdict = build_momentdict(sdp, renamemoments)

    Build a dict `momentdict::SortedDict{Exponent, String}` that stores short names for moments (Exponent) of `sdp`.
    NOTE: should be extended to Variable and (String, Variable) later...
"""
function build_momentdict(sdp, renamemoments::Bool)
    momentdict = SortedDict{Exponent, String}()

    n_moment=0
    n_matvar=0
    n_scalvar=0

    momentdict[Exponent()] = "1"
    ## Treating blocks
    for (moment, blockname, γ, δ) in keys(sdp.blocks)
        α, β = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, α) || (n_moment+=1; momentdict[α] = shortname_moment(n_moment))
            haskey(momentdict, β) || (n_moment+=1; momentdict[β] = shortname_moment(n_moment))
            haskey(momentdict, γ) || (n_moment+=1; momentdict[γ] = shortname_moment(n_moment))
            haskey(momentdict, δ) || (n_moment+=1; momentdict[δ] = shortname_moment(n_moment))
        else
            haskey(momentdict, α) || (momentdict[α] = format_string(α))
            haskey(momentdict, β) || (momentdict[β] = format_string(β))
            haskey(momentdict, γ) || (momentdict[γ] = format_string(γ))
            haskey(momentdict, δ) || (momentdict[δ] = format_string(δ))
        end
    end

    for (moment, var) in keys(sdp.lin)
        α, β = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, α) || (n_moment+=1; momentdict[α] = shortname_moment(n_moment))
            haskey(momentdict, β) || (n_moment+=1; momentdict[β] = shortname_moment(n_moment))
        else
            momentdict[α] = format_string(α)
            momentdict[β] = format_string(β)
        end
    end

    for (moment, blockname, var) in keys(sdp.linsym)
        α, β = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, α) || (n_moment+=1; momentdict[α] = shortname_moment(n_moment))
            haskey(momentdict, β) || (n_moment+=1; momentdict[β] = shortname_moment(n_moment))
        else
            momentdict[α] = format_string(α)
            momentdict[β] = format_string(β)
        end
    end

    for moment in keys(sdp.cst)
        α, β = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, α) || (n_moment+=1; momentdict[α] = shortname_moment(n_moment))
            haskey(momentdict, β) || (n_moment+=1; momentdict[β] = shortname_moment(n_moment))
        else
            momentdict[α] = format_string(α)
            momentdict[β] = format_string(β)
        end
    end

    return momentdict
end

function build_ctrkeysset(sdp)
    ctr_keys = SortedSet{Moment}()

    for ((moment, blockname, γ, δ), λ) in sdp.blocks
        push!(ctr_keys, moment)
    end
    for ((moment, var), λ) in sdp.lin
        push!(ctr_keys, moment)
    end
    for ((moment, blockname, var), λ) in sdp.linsym
        push!(ctr_keys, moment)
    end
    for (moment, λ) in sdp.cst
        push!(ctr_keys, moment)
    end
    return ctr_keys
end


function print_blocksfile(io::IO, sdpblocks::SDPBlocks, momentdict; indentedprint=false)
    println(io, "## Description of the matrices A_ji for the problem:")
    println(io, "##         max     ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0")
    println(io, "##         s.t.    ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  ==  0")
    println(io, "## Constraints keys are j → (j_conj, j_expl, clique).")
    println(io, "## Objective keys are 0 → (1,1, *) for any *.")
    println(io, "#")

    cstrlenα = maximum(x->length(momentdict[x[1].conj_part]), keys(sdpblocks))
    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ = maximum(x->length(momentdict[x[1].expl_part]), keys(sdpblocks))
    cstrlenβ= max(cstrlenβ, length("j_expl"))
    cliquelen = maximum(x->length(x[1].clique), keys(sdpblocks))
    cliquelen= max(cliquelen, length("clique"))
    blocklen = maximum(x->length(x[2]), keys(sdpblocks))
    blocklen= max(blocklen, length("Zi"))
    rowlen = maximum(x->length(momentdict[x[3]]), keys(sdpblocks))
    rowlen = max(rowlen, length("k"))
    collen = maximum(x->length(momentdict[x[4]]), keys(sdpblocks))
    collen = max(collen, length("l"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "Zi", blocklen, indentedprint=indentedprint)
    print_string(io, "k", rowlen, indentedprint=indentedprint)
    print_string(io, "l", collen, indentedprint=indentedprint)
    print_string(io, "Real(A_ij[k, l])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(A_ij[k, l])", 23, indentedprint=indentedprint)
    println(io)

    for ((moment, blockname, γ, δ), λ) in sdpblocks
        α, β = moment.conj_part, moment.expl_part
        print_string(io, momentdict[α], cstrlenα, indentedprint=indentedprint)
        print_string(io, momentdict[β], cstrlenβ, indentedprint=indentedprint)
        print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
        print_string(io, blockname, blocklen, indentedprint=indentedprint)
        print_string(io, momentdict[γ], rowlen, indentedprint=indentedprint)
        print_string(io, momentdict[δ], collen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end


function print_linfile(io::IO, sdplin::SDPLin, sdplinsym::SDPLinSym, momentdict; indentedprint=false)
    println(io, "## Description of the vectors b_j for the problem:")
    println(io, "##         max     ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0")
    println(io, "##         s.t.    ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  ==  0")
    println(io, "## Constraints keys are j → (j_conj, j_expl, clique).")
    println(io, "## Objective keys are 0 → (1,1, *) for any *.")
    println(io, "#")

    cstrlenα = 0
    cstrlenβ = 0
    cliquelen = 0
    varlen = varlensym = 0

    if length(union(keys(sdplin), keys(sdplinsym))) > 0
        cstrlenα = maximum(x->length(momentdict[x[1].conj_part]), union(keys(sdplin), keys(sdplinsym)))
        cstrlenβ = maximum(x->length(momentdict[x[1].expl_part]), union(keys(sdplin), keys(sdplinsym)))
        cliquelen = maximum(x->length(x[1].clique), union(keys(sdplin), keys(sdplinsym)))
        varlen = length(sdplin)!=0 ? maximum(x->length(format_string(x[2])), keys(sdplin)) : 0
        varlensym = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[3], x[2])), keys(sdplinsym)) : 0
    end

    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ= max(cstrlenβ, length("j_expl"))
    cliquelen= max(cliquelen, length("clique"))
    varlen = max(varlen, varlensym, length("x[k]"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "x[k]", varlen, indentedprint=indentedprint)
    print_string(io, "Real(b_j[k])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(b_j[k])", 23, indentedprint=indentedprint)
    println(io)

    if length(sdplin)!=0
        for ((moment, var), λ) in sdplin
            α, β = moment.conj_part, moment.expl_part
            print_string(io, momentdict[α], cstrlenα, indentedprint=indentedprint)
            print_string(io, momentdict[β], cstrlenβ, indentedprint=indentedprint)
            print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
            print_string(io, format_string(var), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
        end
    end
    if length(sdplinsym) != 0
        for ((moment, blockname, var), λ) in sdplinsym
            α, β = moment.conj_part, moment.expl_part
            print_string(io, momentdict[α], cstrlenα, indentedprint=indentedprint)
            print_string(io, momentdict[β], cstrlenβ, indentedprint=indentedprint)
            print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
            print_string(io, format_string(var, blockname), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
        end
    end
end


function print_cstfile(io::IO, sdpcst::SDPCst, momentdict, ctr_keys::SortedSet{Moment}; indentedprint=false)
    println(io, "## Description of the scalars c_j for the problem:")
    println(io, "##         max     ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0")
    println(io, "##         s.t.    ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  ==  0")
    println(io, "## Constraints keys are j → (j_conj, j_expl, clique).")
    println(io, "## Objective keys are 0 → (1,1, *) for any *.")
    println(io, "#")

    cstrlenα = maximum(x->length(momentdict[x.conj_part]), ctr_keys)
    cstrlenα= max(cstrlenα, length("#j_conj"))
    cstrlenβ = maximum(x->length(momentdict[x.expl_part]), ctr_keys)
    cstrlenβ= max(cstrlenβ, length("j_expl"))
    cliquelen = maximum(x->length(x.clique), ctr_keys)
    cliquelen= max(cliquelen, length("clique"))

    print_string(io, "#j_conj", cstrlenα, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlenβ, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "Real(c_j)", 23, indentedprint=indentedprint)
    print_string(io, "Imag(c_j)", 23, indentedprint=indentedprint)
    println(io)

    for moment in ctr_keys
        α, β = moment.conj_part, moment.expl_part
        λ = haskey(sdpcst, moment)?sdpcst[moment]:0
        print_string(io, momentdict[α], cstrlenα, indentedprint=indentedprint)
        print_string(io, momentdict[β], cstrlenβ, indentedprint=indentedprint)
        print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(λ), imag(λ))
    end
end


function print_typesfile(io::IO, block_to_vartype)
    println(io, "## Description of the matrix and scalar variables Zi and x[k] for the problem:")
    println(io, "##         max     ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0")
    println(io, "##         s.t.    ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  ==  0")
    println(io, "## Matrix variable types Zi are \"SDPC\" or \"SDP\".")
    println(io, "## By default, scalar variables x[k] are assumed to be real, free.")
    println(io, "#")

    cstrlen = maximum(x->length(x), keys(block_to_vartype))
    cstrlen = max(cstrlen, length("#Zi"))

    print_string(io, "#Zi", cstrlen, alignright=false); println(io, " type")

    for (blockname, vartype) in block_to_vartype
        if vartype in Set([:SDP, :SDPC])
            print_string(io, blockname, cstrlen)
            println(io, " $(string(vartype))")
        end
        # warn("Ignoring variable $blockname of type $vartype") ## NOTE: better logging system.
    end
end


function print_namesfile(io::IO, momentdict)
    println(io, "## Description of the scalars c_j for the problem:")
    println(io, "##         max     ∑ A_0i[k,l] × Zi[k,l] + ∑ b_0[k] × x[k] + c_0")
    println(io, "##         s.t.    ∑ A_ji[k,l] × Zi[k,l] + ∑ b_j[k] × x[k] + c_j  ==  0")
    println(io, "## Constraints keys are j → (j_conj, j_expl, clique).")
    println(io, "## Objective keys are 0 → (1,1, *) for any *.")
    println(io, "#")
    println(io, "#shortname  Explicit_name")

    for (α, shortname) in momentdict
        println(io, "$shortname $(format_string(α))")
    end
end