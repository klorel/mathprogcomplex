"""
    m, Zi, yα_re, yα_im = make_JuMPproblem(SDP_SOS::SDPPrimal, mysolver)

    Convert the SDP_SOS problem into a `m` JuMP problem, with primal variables `Zi` and dual variables `Zi`
"""
function make_JuMPproblem(SDP::SDPPrimal, mysolver)
    m = Model(solver = mysolver)

    
    ## Variables
    Zi = Dict{String, Array{JuMP.Variable, 2}}()
    for (cstrname, n) in SDP_SOS.variables
        var = @variable(m, [1:2*n,1:2*n], SDP, basename=cstrname)
        Zi[cstrname] = var
        cstr1 = @constraint(m, [i=1:n,j=1:n], var[i, j] == var[i+n, j+n])
        for i in 1:n, j in 1:n
            println("$i, $j  --> $(cstr1[i, j])")
        end
        cstr2 = @constraint(m, [i=1:n,j=1:i], var[n+i, j] == - var[n+j, i])
        for i in 1:n, j in 1:i
            println("$i, $j --> $(cstr2[i, j])")
        end
    end

    
    ## Setting objective
    expo = Exponent()

    objectif_re, _ = formconstraint(SDP_SOS.objective, Zi)
    @objective(m, Max, objectif_re)


    ## Setting constraints
    d = relax_ctx.di["moment_cstr"]
    vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)

    #Constraint storage
    @constraintref yα_re[1:length(realexpos), 1:length(realexpos)]
    @constraintref yα_im[1:length(conjexpos), 1:length(conjexpos)]

    expo2int = Dict{Exponent, Int}()
    i = 1
    for expo in sort(collect(realexpos))
        expo2int[expo] = i
        i += 1
    end

    yα_re[1,1] = @constraint(m, 1==1)
    yα_im[1,1] = @constraint(m, 1==1)

    #Building constraints
    for re in realexpos, im in conjexpos
        expo = product(re, im)

        ## NOTE: (alpha, beta) and (beta, alpha) should provide the same constraint, to be checked.
        if expo != Exponent()
            println("--> $(expo2int[re]), $(expo2int[conj(im)])")
            l_re, l_im = formconstraint(SDP_SOS.constraints[expo], Zi)

            if l_re != 0
                yα_re[expo2int[re], expo2int[conj(im)]] = @constraint(m, l_re == 0)
                println(yα_re[expo2int[re], expo2int[conj(im)]])
            else
                println("-> Expo $expo provided no constraint")
            end

            if l_im != 0
                yα_im[expo2int[re], expo2int[conj(im)]] = @constraint(m, l_im == 0)
                println(yα_im[expo2int[re], expo2int[conj(im)]])
            else
                println("-> Expo $expo provided no constraint")
            end
        end
    end

    return m, Zi, yα_re, yα_im

end


function Λreal(A::AbstractMatrix)
    n, m = size(A)
    return A[1:Int(n/2), 1:Int(m/2)]
end

function Λimag(A::AbstractMatrix)
    n, m = size(A)
    return A[(Int(n/2+1)):n, 1:Int(m/2)]
end


function formconstraint(SDPform, Zi)
    l_re = real(SDPform.scal)
    l_im = imag(SDPform.scal)
    for (i, Biα) in SDPform.mats
        l_re -= vecdot(Λreal(Zi[i]), real(Biα)) - vecdot(Λimag(Zi[i]), imag(Biα))
        l_im -= vecdot(Λreal(Zi[i]), imag(Biα)) + vecdot(Λimag(Zi[i]), real(Biα))
    end
    return l_re, l_im
end