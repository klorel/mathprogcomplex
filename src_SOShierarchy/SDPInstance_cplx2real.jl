
function SDPInstance_cplx2real(sdp::SDPInstance)
    sdpblocks = SDPBlocks()
    sdplinsym = SDPLinSym()
    sdplin = SDPLin()
    sdpcst = SDPCst()
    block_to_vartype = SortedDict{String, Symbol}()

    ## Complex blocks to real
    for ((moment, block_name, γ, δ), coeff) in sdp.blocks
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        γ_re, γ_im = cplx2real_sdpctr(γ)
        δ_re, δ_im = cplx2real_sdpctr(δ)


        println("$moment  is $ctr_re, $ctr_im")
        # Convert complex linear term to real ones
        haskey(sdpblocks, (ctr_re, block_name, γ_re, δ_re)) && warn("haskey(sdpblocks, ($ctr_re, $block_name, $γ_re, $δ_re))")
        sdpblocks[(ctr_re, block_name, γ_re, δ_re)] = real(coeff)
        info("adding sdpblocks -> ($ctr_re, $block_name, $γ_re, $δ_re))")

        haskey(sdpblocks, (ctr_re, block_name, γ_im, δ_im)) && warn("haskey(sdpblocks, ($ctr_re, $block_name, $γ_im, $δ_im))")
        sdpblocks[(ctr_re, block_name, γ_re, δ_im)] = -imag(coeff)
        info("adding sdpblocks -> ($ctr_re, $block_name, $γ_im, $δ_im))")

        haskey(sdpblocks, (ctr_im, block_name, γ_re, δ_re)) && warn("haskey(sdpblocks, ($ctr_im, $block_name, $γ_re, $δ_re))")
        sdpblocks[(ctr_im, block_name, γ_re, δ_re)] = imag(coeff)
        info("adding sdpblocks -> ($ctr_im, $block_name, $γ_re, $δ_re))")

        haskey(sdpblocks, (ctr_im, block_name, γ_im, δ_im)) && warn("haskey(sdpblocks, ($ctr_im, $block_name, $γ_im, $δ_im))")
        sdpblocks[(ctr_im, block_name, γ_re, δ_im)] = real(coeff)
        info("adding sdpblocks -> ($ctr_im, $block_name, $γ_im, $δ_im))")

        # Add constraint on matrix variable: real part symmetric & 1,1 == 2,2
        Xi_ctrname_re = get_Xictrname_re(block_name, γ, δ)

        haskey(sdpblocks, (Xi_ctrname_re, block_name, γ_re, δ_re)) && warn("hasekey sdpblocks[($(Xi_ctrname_re), $block_name, $γ_re, $δ_re)]")
        sdpblocks[(Xi_ctrname_re, block_name, γ_re, δ_re)] = 1    #haskey(sdpblocks, (Xi_ctrname_re, block_name, γ_re, δ_re)) ||
        print_with_color(:light_cyan, "adding sdpblocks -> ($(Xi_ctrname_re), $block_name, $γ_re, $δ_re)\n")

        haskey(sdpblocks, (Xi_ctrname_re, block_name, γ_im, δ_im)) && warn("hasekey sdpblocks[($(Xi_ctrname_re), $block_name, $γ_im, $δ_im)]")
        sdpblocks[(Xi_ctrname_re, block_name, γ_im, δ_im)] = -1    #haskey(sdpblocks, (Xi_ctrname_re, block_name, γ_im, δ_im)) ||
        print_with_color(:light_cyan, "adding sdpblocks -> ($(Xi_ctrname_re), $block_name, $γ_im, $δ_im)\n")

        # Add constraint on matrix variable: imad part antisymmetric && 1,2 == 2,1^T
        Xi_ctrname_im = get_Xictrname_im(block_name, γ, δ)

        haskey(sdpblocks, (Xi_ctrname_im, block_name, γ_re, δ_im)) && warn("hasekey sdpblocks[($(Xi_ctrname_im), $block_name, $γ_re, $δ_im)]")
        sdpblocks[(Xi_ctrname_im, block_name, γ_re, δ_im)] = 1
        print_with_color(:light_cyan, "adding sdpblocks -> ($(Xi_ctrname_im), $block_name, $γ_re, $δ_im)\n")

        haskey(sdpblocks, (Xi_ctrname_im, block_name, γ_im, δ_re)) && warn("hasekey sdpblocks[($(Xi_ctrname_im), $block_name, $γ_im, $δ_re)]")
        sdpblocks[(Xi_ctrname_im, block_name, γ_im, δ_re)] = 1
        print_with_color(:light_cyan, "adding sdpblocks -> ($(Xi_ctrname_im), $block_name, $γ_im, $δ_re)\n")
    end

    ## Complex symetric blocks to real
    for ((moment, block_name, var), coeff) in sdp.linsym
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        var_re, var_im = cplx2real_sdpctr(var)

        sdplinsym[(ctr_re, block_name, var_re)] = real(coeff)
        sdplinsym[(ctr_re, block_name, var_im)] = -imag(coeff)

        sdplinsym[(ctr_im, block_name, var_re)] = imag(coeff)
        sdplinsym[(ctr_im, block_name, var_im)] = real(coeff)
    end

    ## Complex vectors to real
    for ((moment, var), coeff) in sdp.lin
        ctr_re, ctr_im = cplx2real_sdpctr(moment)
        var_re, var_im = cplx2real_sdpctr(var)

        sdplin[(ctr_re, var_re)] = real(coeff)
        sdplin[(ctr_re, var_im)] = -imag(coeff)

        sdplin[(ctr_im, var_re)] = imag(coeff)
        sdplin[(ctr_im, var_im)] = real(coeff)
    end

    ## Complex coeff to real
    for (moment, coeff) in sdp.cst
        ctr_re, ctr_im = cplx2real_sdpctr(moment)

        sdpcst[ctr_re] = real(coeff)
        sdpcst[ctr_im] = imag(coeff)
    end

    ## Convert variable types
    for (block, vartype) in sdp.block_to_vartype
        if vartype == :SDPC
            block_to_vartype[block] = :SDP
        elseif vartype == :SymC
            block_to_vartype[block] = :Sym
        else
            error("SDPInstance_cplx2real(): Unhandled matrix type $vartype for $block")
        end
    end

    return SDPInstance(block_to_vartype, sdpblocks, sdplinsym, sdplin, sdpcst)
end



function cplx2real_sdpctr(moment::Moment)
    ctr_re = Moment(moment.conj_part, moment.expl_part, moment.clique*"_Re")
    ctr_im = Moment(moment.conj_part, moment.expl_part, moment.clique*"_Im")

    return ctr_re, ctr_im
end

# function cplx2real_sdpctr(expo::Exponent)
#     expo_re, expo_im = Exponent(), Exponent()
#     expo_re = Exponent(Variable(string(expo, "_Re"), Real))
#     expo_im = Exponent(Variable(string(expo, "_Im"), Real))
#     return expo_re, expo_im
# end

function get_Xictrname_re(block_name::String, γ::Exponent, δ::Exponent)
    return Moment(γ, δ, block_name*"_ReCtr")
    # return Exponent(Variable(block_name*"_Re", Real)), Exponent(Variable(string(γ, "_", δ), Real))
end

function get_Xictrname_im(block_name::String, γ::Exponent, δ::Exponent)
    return Moment(γ, δ, block_name*"_ImCtr")
    # return Exponent(Variable(block_name*"_Im", Real)), Exponent(Variable(string(γ, "_", δ), Real))
end