function export_SDP(relax_ctx, sdpbody, sdprhs, path)
    if relax_ctx.hierarchykind != :Real
        error("export_SDP(): Only real hierarchy for now... C -> R conversion still to be done")
    end

    
    # Build body and rhs proper bodys
    body = SortedDict{Tuple{String, String, String, String}, Float64}()
    rhs = SortedDict{String, Float64}()

    for ((cstrname, blocname, α, β), Bi) in sdpbody
        cstrstr = format_string(α, β)
        blocstr = format_string(cstrname, blocname)
        for ((γ, δ), λ) in Bi
            i = format_string(γ)
            j = format_string(δ)
            if j <= i
                body[(cstrstr, blocstr, i, j)] = λ
            end
        end
    end

    for ((α, β), fαβ) in sdprhs
        cstrstr = format_string(α, β)
        rhs[cstrstr] = fαβ
    end

    # Export blocks of constraints
    fbody = open(joinpath(path, "body.sdp"), "w")

    cstrlen = maximum(x->length(x[1]), keys(body))
    cstrlen = max(cstrlen, length("# cstrname"))
    bloclen = maximum(x->length(x[2]), keys(body))
    bloclen = max(bloclen, length("blocname"))
    expo1len = maximum(x->length(x[3]), keys(body))
    expo1len = max(expo1len, length("row"))
    expo2len = maximum(x->length(x[3]), keys(body))
    expo2len = max(expo2len, length("col"))

    print_string(fbody, "# cstrname", cstrlen)
    print_string(fbody, "blocname", bloclen)
    print_string(fbody, "row", expo1len)
    print_string(fbody, "col", expo2len)
    println(fbody, " val")

    for ((cstrname, blocname, i, j), val) in body
        print_string(fbody, cstrname, cstrlen)
        print_string(fbody, blocname, bloclen)
        print_string(fbody, i, expo1len)
        print_string(fbody, j, expo2len)
        println(fbody, " $val")
    end
    close(fbody)

    # TODO: Export linear terms of constraints
    # cstr -> var -> val

    # Export rhs
    frhs = open("rhs.sdp", "w")
    cstrlen = maximum(x->length(x), keys(rhs))

    cstrlen = max(cstrlen, length("cstrname")+1)
    print(frhs, "#")
    print_string(frhs, "cstrname", cstrlen-1)
    println(frhs, " rhs_val")
    for (cstrname, val) in rhs
        print_string(frhs, cstrname, cstrlen)
        println(frhs, " $val")
    end
    close(frhs)

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