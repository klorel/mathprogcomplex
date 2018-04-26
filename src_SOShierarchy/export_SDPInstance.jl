function export_SDP(relax_ctx, sdp::SDPInstance, path)
    if relax_ctx.hierarchykind != :Real
        error("export_SDP(): Only real hierarchy for now... C -> R conversion still to be done")
    end

    

    # Export blocks of constraints
    blocks_file = joinpath(path, "blocks.sdp")
    !isfile(blocks_file) || rm(blocks_file)

    fblocks = open(blocks_file, "w")
    print(fblocks, sdp.blocks)
    close(fblocks)

    # Export linear part of constraints
    lin_file = joinpath(path, "lin.sdp")
    !isfile(lin_file) || rm(lin_file)

    if length(sdp.lin) > 0
        flin = open(lin_file, "w")
        print(flin, sdp.lin)
        close(flin)
    end

    # Export constants of constraints
    cst_file = joinpath(path, "const.sdp")
    !isfile(cst_file) || rm(cst_file)

    cst = open(cst_file, "w")
    print(cst, sdp.cst)
    close(cst)


    # Export bloc types
    ftypes = open("types.sdp", "w")
    cstrlen = maximum(x->length(x), keys(relax_ctx.cstrtypes))
    cstrlen = max(cstrlen, length("# cstrname"))
    print_string(ftypes, "# cstrname", cstrlen)
    println(ftypes, " var_type")
    for (cstrname, cstrtype) in relax_ctx.cstrtypes
        print_string(ftypes, cstrname, cstrlen)
        println(ftypes, " $(string(cstrtype))")
    end
    close(ftypes)
end