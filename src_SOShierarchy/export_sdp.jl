function export_SDP(relax_ctx, sdp::SDPInstance, path)
    if relax_ctx.hierarchykind != :Real
        error("export_SDP(): Only real hierarchy for now... C -> R conversion still to be done")
    end

    
    # Build blocks and rhs proper
    blocks = SortedDict{Tuple{String, String, String, String}, Float64}()
    
    for ((cstrname, blocname, α, β), Bi) in sdp.blocks
        cstrstr = format_string(α, β)
        blocstr = format_string(cstrname, blocname)
        for ((γ, δ), λ) in Bi
            i = format_string(γ)
            j = format_string(δ)
            # if j <= i
            @assert typeof(λ) == Float64
            @assert !haskey(blocks, (cstrstr, blocstr, i, j))
            blocks[(cstrstr, blocstr, i, j)] = λ
        end
    end

    lin = SortedDict{Tuple{String, String}, Float64}()
    for ((α, β, α), λ) in sdp.lin
        cstrstr = format_string(α, β)
        i = format_string(α)
        @assert typeof(λ) == Float64
        @assert !haskey(lin, (cstrstr, i))
        lin[(cstrstr, i)] = λ
    end

    cnst = SortedDict{String, Float64}()
    for ((α, β), fαβ) in sdp.cnst
        cstrstr = format_string(α, β)
        cnst[cstrstr] = fαβ
    end

    # Export blocks of constraints
    fblocks = open(joinpath(path, "blocks.sdp"), "w")

    cstrlen = maximum(x->length(x[1]), keys(blocks))
    cstrlen = max(cstrlen, length("# cstrname"))
    bloclen = maximum(x->length(x[2]), keys(blocks))
    bloclen = max(bloclen, length("blocname"))
    expo1len = maximum(x->length(x[3]), keys(blocks))
    expo1len = max(expo1len, length("row"))
    expo2len = maximum(x->length(x[3]), keys(blocks))
    expo2len = max(expo2len, length("col"))

    print_string(fblocks, "# cstrname", cstrlen)
    print_string(fblocks, "blocname", bloclen)
    print_string(fblocks, "row", expo1len)
    print_string(fblocks, "col", expo2len)
    println(fblocks, " val")

    for ((cstrname, blocname, i, j), val) in blocks
        print_string(fblocks, cstrname, cstrlen)
        print_string(fblocks, blocname, bloclen)
        print_string(fblocks, i, expo1len)
        print_string(fblocks, j, expo2len)
        println(fblocks, " $val")
    end
    close(fblocks)

    # Export linear
    if length(lin) > 0
        flin = open("lin.sdp", "w")

        cstrlen = maximum(x->length(x[1]), keys(lin))
        cstrlen = max(cstrlen, length("cstrname")+1)
        varlen = maximum(x->length(x[2]), keys(lin))
        varlen = max(varlen, length("varname"))

        print_string(flin, "# cstrname", cstrlen)
        print_string(flin, "varlen", varlen)
        println(flin, " val")
        
        for ((cstrname, varname), val) in lin
            print_string(flin, cstrname, cstrlen)
            print_string(flin, varname, varlen)
            println(flin, " $val")
        end
        close(flin)
    else
        touch("lin.sdp")
    end

    # Export const
    fconst = open("const.sdp", "w")
    
    cstrlen = maximum(x->length(x), keys(cnst))
    cstrlen = max(cstrlen, length("# cstrname"))

    print_string(fconst, "# cstrname", cstrlen)
    println(fconst, " const_val")
    for (cstrname, val) in cnst
        print_string(fconst, cstrname, cstrlen)
        println(fconst, " $val")
    end
    close(fconst)

    # Export bloc types
    ftypes = open("types.sdp", "w")
    cstrlen = maximum(x->length(x), keys(relax_ctx.cstrtypes))
    print(ftypes, "#")
    print_string(ftypes, "cstrname", cstrlen-1)
    println(ftypes, " var_type")
    for (cstrname, cstrtype) in relax_ctx.cstrtypes
        print_string(ftypes, cstrname, cstrlen)
        println(ftypes, " $(string(cstrtype))")
    end
    close(ftypes)
end


function format_string(α::Exponent, β::Exponent)
    s = "$α,$β"
    return replace(s, " ", "_")
end

function format_string(s1::String, s2::String)
    s = "$s1,$s2"
    return replace(s, " ", "_")
end

function format_string(α::Exponent)
    s = "$α"
    return replace(s, " ", "_")
end