"""
    hassym = has_phasesymmetry(relax.ctx, pb)

    Check that the objective and all constraints present a phase-shift (complex torus) invariance.
"""
function has_phasesymmetry(relax_ctx::RelaxationContext, pb::Problem)
    polykind = relax_ctx.hierarchykind
    hassym = is_homogeneous(pb.objective, polykind)
    println("$hassym  \t $(pb.objective)")
    for (cstrname, cstr) in pb.constraints
        hassym = hassym && is_homogeneous(cstr.p, polykind)
        println("$hassym  \t $(cstr.p)")
    end
    return hassym
end


"""
enforce_phaseinvariance!(relax_ctx, momentmatrices)
"""
function enforce_phaseinvariance!(relax_ctx::RelaxationContext, momentmatrices::Dict{Tuple{String, String}, MomentMatrix})
    # Generate appropriate point from the moment contraint
    vars = momentmatrices[("moment_cstr", "oneclique")].vars
    d = momentmatrices[("moment_cstr", "oneclique")].order

    println("enforce_phaseinvariance!(): pt_null construction...")

    realexpos = compute_exponents(vars, d)
    conjexpos = compute_exponents(vars, d, compute_conj=true)
    pt_null = Point()
    for realexpo in realexpos, conjexpo in conjexpos
        expo = product(realexpo, conjexpo)
        println(expo)
        if !is_homogeneous(expo, relax_ctx.hierarchykind)
            pt_null.coords[expo] = 0
        end
    end

    println("pt_null: $pt_null")

    for (key, mmt) in momentmatrices
        # evaluate the moment matrix at that point
        println("Dealing with $key...")
        println("Before: $mmt")
        momentmatrices[key] = evaluate(mmt, pt_null) # TODO : evaluate!
        println("After: $(momentmatrices[key])")
    end
end