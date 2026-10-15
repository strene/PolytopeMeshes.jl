# Conditions on background sites that guarantee conformity of the Voronoi
# diagram to face and cell constraints.

"""
    surface_sufficient_condition(pts, centers, radii; tol)

Remove the points strictly inside any of the circles/balls given by `centers`
and `radii` (UPR: `surfaceSufCond2D`/`surfaceSufCond3D`). Returns
`(kept_points, removed)` where `removed` is a boolean mask.
"""
function surface_sufficient_condition(pts::AbstractVector{SVector{D, Float64}}, centers::AbstractVector{SVector{D, Float64}},
        radii::AbstractVector{<:Real}; tol = nothing) where {D}
    removed = falses(length(pts))
    if isempty(pts) || isempty(centers)
        return collect(pts), removed
    end
    if isnothing(tol)
        lo, hi = bounding_box(centers)
        tol = 10 * eps() * max(norm(hi - lo), maximum(radii))
    end
    tree = KDTree(collect(pts))
    for (c, R) in zip(centers, radii)
        R2 = R^2 - tol
        R2 <= 0 && continue
        for i in inrange(tree, c, R)
            if sum(abs2, pts[i] - c) < R2
                removed[i] = true
            end
        end
    end
    return pts[.!removed], removed
end

"""
    remove_conflict_points(P1, P2, dist)

Remove the points in `P1` closer than `dist[j]` to `P2[j]` for any `j` (UPR:
`removeConflictPoints`). Returns `(kept_points, removed)`.
"""
function remove_conflict_points(P1::AbstractVector{SVector{D, Float64}}, P2::AbstractVector{SVector{D, Float64}}, dist) where {D}
    removed = falses(length(P1))
    if isempty(P1) || isempty(P2)
        return collect(P1), removed
    end
    dist = dist isa Real ? fill(Float64(dist), length(P2)) : dist
    tree = KDTree(collect(P1))
    for (p, d) in zip(P2, dist)
        for i in inrange(tree, p, d)
            if norm(P1[i] - p) < d
                removed[i] = true
            end
        end
    end
    return P1[.!removed], removed
end

"""
    line_sufficient_condition(pts, lines; tol)

Remove points inside the balls having consecutive sites of a cell constraint
as diameter, so that consecutive constraint sites stay Voronoi neighbors (UPR:
`lineSufCond3D`). `lines` is a vector of vectors of consecutive sites.
Returns `(kept_points, removed)`.
"""
function line_sufficient_condition(pts::AbstractVector{SVector{D, Float64}}, lines; tol = 10 * eps()) where {D}
    centers = SVector{D, Float64}[]
    radii = Float64[]
    for l in lines
        for k in 1:(length(l) - 1)
            push!(centers, (l[k] + l[k + 1]) / 2)
            push!(radii, norm(l[k + 1] - l[k]) / 2)
        end
    end
    removed = falses(length(pts))
    if isempty(pts) || isempty(centers)
        return collect(pts), removed
    end
    tree = KDTree(collect(pts))
    for (c, R) in zip(centers, radii)
        for i in inrange(tree, c, R + sqrt(tol))
            if sum(abs2, pts[i] - c) < R^2 + tol
                removed[i] = true
            end
        end
    end
    return pts[.!removed], removed
end
