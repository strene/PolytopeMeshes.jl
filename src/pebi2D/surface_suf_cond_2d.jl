"""
Enforce the sufficient surface condition.
Ported from MRST's surfaceSufCond2D.
"""

"""
    surface_suf_cond_2d(p, F)

Remove points that violate the sufficient surface condition, i.e., points
that lie inside any of the circles used to create the fault sites.

# Arguments
- `p`: n×2 array of candidate points
- `F`: SurfaceSitesResult from `surface_sites_2d`

# Returns
- `p_filtered`: Points that satisfy the condition
- `removed`: Boolean array indicating removed points
"""
function surface_suf_cond_2d(p::Matrix{Float64}, F::SurfaceSitesResult)
    if size(F.c_CC, 1) == 0 || size(p, 1) == 0
        return copy(p), falses(size(p, 1))
    end

    TOL = 10 * (maximum(F.f_pts) - minimum(F.f_pts)) * eps()
    nc = size(F.c_CC, 1)
    np = size(p, 1)

    removed = falses(np)
    for i in 1:nc
        cr_sqr = F.c_R[i]^2
        for j in 1:np
            if !removed[j]
                dist_sqr = (F.c_CC[i, 1] - p[j, 1])^2 + (F.c_CC[i, 2] - p[j, 2])^2
                if dist_sqr < cr_sqr - TOL
                    removed[j] = true
                end
            end
        end
    end

    return p[.!removed, :], removed
end
