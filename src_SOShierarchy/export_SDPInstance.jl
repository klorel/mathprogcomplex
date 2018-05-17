function export_SDP(relax_ctx, sdp::SDPInstance, path)
    if relax_ctx.hierarchykind != :Real
        error("export_SDP(): Only real hierarchy for now... C -> R conversion still to be done")
    end

    # Export blocks of constraints
    blocks_file = joinpath(path, "blocks.sdp")
    !isfile(blocks_file) || rm(blocks_file)

    fblocks = open(blocks_file, "a")
    println(fblocks, "## Description of the matrices A_ji for the problem:")
    println(fblocks, "##         max     ∑ A_0i ⋅ Zi + b_0 ⋅ x + c0")
    println(fblocks, "##         s.t.    ∑ A_ji ⋅ Zi + b_j ⋅ x + c_j  ==  0")
    println(fblocks, "## Objective key is \"1,1\"")
    println(fblocks, "#")
    print(fblocks, sdp.blocks)
    close(fblocks)

    # Export linear part of constraints
    lin_file = joinpath(path, "lin.sdp")
    !isfile(lin_file) || rm(lin_file)

    flin = open(lin_file, "a")
    println(flin, "## Description of the vectors b_j for the problem:")
    println(flin, "##         max     ∑ A_0i ⋅ Zi + b_0 ⋅ x + c0")
    println(flin, "##         s.t.    ∑ A_ji ⋅ Zi + b_j ⋅ x + c_j  ==  0")
    println(flin, "## Objective key is \"1,1\"")
    println(flin, "#")
    print(flin, sdp.lin)
    close(flin)

    # Export constants of constraints
    cst_file = joinpath(path, "const.sdp")
    !isfile(cst_file) || rm(cst_file)

    fcst = open(cst_file, "a")
    println(fcst, "## Description of the scalars c_j for the problem:")
    println(fcst, "##         max     ∑ A_0i ⋅ Zi + b_0 ⋅ x + c0")
    println(fcst, "##         s.t.    ∑ A_ji ⋅ Zi + b_j ⋅ x + c_j  ==  0")
    println(fcst, "## Objective key is \"1,1\"")
    println(fcst, "#")
    print(fcst, sdp.cst)
    close(fcst)


    # Export bloc types
    types_file = joinpath(path, "types.sdp")
    !isfile(types_file) || rm(types_file)

    ftypes = open(types_file, "a")
    println(ftypes, "## Description of the matrix variables Zi for the problem:")
    println(ftypes, "##         max     ∑ A_0i ⋅ Zi + b_0 ⋅ x + c0")
    println(ftypes, "##         s.t.    ∑ A_ji ⋅ Zi + b_j ⋅ x + c_j  ==  0")
    println(ftypes, "## Matrix type are \"SDP\" or \"Sym\".")
    println(ftypes, "#")
    cstrlen = maximum(x->length(x), keys(sdp.block_to_vartype))
    cstrlen = max(cstrlen, length("# Matrix variable key i"))
    print_string(ftypes, "# Matrix variable key i", cstrlen); println(ftypes, " # Matrix type")
    for (blockname, vartype) in sdp.block_to_vartype
        print_string(ftypes, blockname, cstrlen)
        println(ftypes, " $(string(vartype))")
    end
    close(ftypes)
end