"""
    mm = MomentMatrix(vars::Set{Variable}, d)

    Build the moment matrix corresponding to the moment of degree up to `d` of the `vars` polynomial algebra.
"""
function MomentMatrix(vars::Set{Variable}, d::Int)
    mm = Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    for cexp in conjexpos
        for rexp in realexpos
            mm[(cexp, rexp)] = cexp*rexp
        end
    end
    return MomentMatrix(mm, copy(vars), d)
end

function copy(mm::MomentMatrix)
    return MomentMatrix(copy(mm.mm), mm.vars, mm.order)
end

function print(io::IO, mm::MomentMatrix)
    for (key, val) in mm.mm
        println(io, "($(key[1]), $(key[2])) ⟶  $val")
    end
end


##########################
# Moment matrix algebra
function product!(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    for (key, val) in mm.mm
        mm.mm[key] = p * mm.mm[key]
    end
end
function product!(λ::Number, mm::MomentMatrix)
    for (key, val) in mm.mm
        mm.mm[key] = λ*val
    end
end

function *(p::T, mm::MomentMatrix) where T<:AbstractPolynomial
    mm_copy = copy(mm)
    product!(p, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, p::T) where T<:AbstractPolynomial = p*mm

function *(λ::Number, mm::MomentMatrix)
    mm_copy = copy(mm)
    product!(λ, mm_copy)
    return mm_copy
end
*(mm::MomentMatrix, λ::Number) = λ*mm

function evaluate(mm::MomentMatrix, pt::Point)
    mm_eval = Dict{Tuple{Exponent, Exponent}, AbstractPolynomial}()
    for (key, p) in mm.mm
        mm_eval[key] = evaluate(p, pt)
    end
    return MomentMatrix(mm_eval, setdiff(mm.vars, keys(pt)), mm.order)
end

#######################################
# Conversion to B base

function MomentMatrixBasis(vars, d, k)
    expo2int = Dict{Exponent, Int}()
    int2expo = Dict{Int, Exponent}()
    i = 1
    for expo in sort(collect(compute_exponents(vars, d-k)))
        expo2int[expo] = i
        int2expo[i] = expo
        i += 1
    end

    return MomentMatrixBasis(Dict{Exponent, AbstractMatrix}(), expo2int, int2expo, length(expo2int))
end


"""
    mmb = convertMMtobase(mm::MomentMatrix, d, k)

    Compute the projection of the `mm` localizing matrix of the `k`-degree constraint on the order `d` moment basis.
    Yield a MomentMatrixBasis object with *dense* matrices.
"""
function convertMMtobase(mm::MomentMatrix, d, k)
    mmb = MomentMatrixBasis(mm.vars, d, k)

    expo2CSCmat = Dict()
    for (key, poly) in mm.mm
        for (expo, λ) in poly
            if (expo.degree.explvar > d) || (expo.degree.conjvar > d)
                warn("convertMMtobase(): Found exponent of degree $(expo.degree) > $d ($expo, at $key of MM matrix)")
            end
            !isnan(λ) || warn("convertMMtobase(): isNaN $key - $expo")

            if !haskey(expo2CSCmat, expo)
                expo2CSCmat[expo] = (Int[], Int[], Complex128[])
            end
            push!(expo2CSCmat[expo][1], mmb.expo2int[conj(key[1])])
            push!(expo2CSCmat[expo][2], mmb.expo2int[key[2]])
            push!(expo2CSCmat[expo][3], λ)
        end
    end

    for (expo, CSCmat) in expo2CSCmat
        mmb.basis[expo] = sparse(CSCmat[1], CSCmat[2], CSCmat[3], mmb.msize, mmb.msize)
    end

    return mmb
end



"""
    momentmatrices = compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)

    Compute the moment and localizing matrices associated with the problem constraints and clique decomposition.
"""
function compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)
    println("\n=== compute_momentmat(problem, max_cliques, cliquevarsbycstr, orderbyclique, relax_ctx)")
    println("Compute the moment and localizing matrices associated with the problem constraints and vlique decomposition.")

    momentmatrices = Dict{String, MomentMatrix}()

    for (cstrname, cstr) in problem.constraints
        vars = Set([Variable(varname, vartype) for (varname, vartype) in problem.variables])
        di, ki = relax_ctx.di[cstrname], relax_ctx.ki[cstrname]
        momentmatrices[cstrname] = MomentMatrix(vars, di - ki) * cstr.p
    end

    println("xxxxx which indicators ? xxxxx")

    return momentmatrices
end
