"""
    hassym = has_phasesymmetry(relax.ctx, pb)

    Check that the objective and all constraints present a phase-shift (complex torus) invariance.
"""
function has_phasesymmetry(relax_ctx::RelaxationContext, pb::Problem)
    polykind = relax_ctx.hierarchykind
    hassym = is_homogeneous(pb.objective, polykind)
    for (cstrname, cstr) in pb.constraints
        hassym = hassym && is_homogeneous(cstr.p, polykind)
    end
    return hassym
end


"""
    enforce_phaseinvariance!(relax_ctx, momentmatrices)

    Remove all the non homogeneous oments from the moment-matrices.
"""
function enforce_phaseinvariance!(relax_ctx::RelaxationContext, momentmatrices::Dict{Tuple{String, String}, MomentMatrix})
    # Remove all the non-homogeneous moments

    for (key, mmt) in momentmatrices
        for ((α, β), p) in mmt.mm
            make_homogeneous!(p, relax_ctx.hierarchykind)
            if p==Polynomial()
                delete!(mmt.mm, (α, β))
            end
        end
    end
end