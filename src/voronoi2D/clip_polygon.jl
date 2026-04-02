"""
Polygon clipping using Sutherland-Hodgman style algorithm.
Ported from MRST's clipPolygon.
"""

"""
    clip_polygon(poly, normals, points)

Clip a convex polygon against a set of half-planes defined by normals and points.
Each half-plane keeps the region where `(x - point) · normal <= 0`.

# Arguments
- `poly`: k×2 matrix of polygon vertices
- `normals`: m×2 matrix of half-plane normals
- `points`: m×2 matrix of points on each half-plane

# Returns
- Clipped polygon vertices (n×2 matrix), or empty matrix if fully clipped
"""
function clip_polygon(poly::Matrix{Float64}, normals::Matrix{Float64},
    points::Matrix{Float64})
    TOL = 100 * eps()

    p = copy(poly)

    for i in 1:size(normals, 1)
        if size(p, 1) < 3
            return zeros(0, 2)
        end

        n = normals[i, :]
        x0 = points[i, :]

        # Calculate signed distance from each vertex to the half-plane
        d = (p .- x0') * n

        if all(d .> TOL)
            return zeros(0, 2)
        elseif all(d .< TOL)
            continue
        end

        # Sutherland-Hodgman clipping
        new_poly = Vector{Vector{Float64}}()
        np = size(p, 1)

        for j in 1:np
            k = mod(j, np) + 1
            d_j = d[j]
            d_k = d[k]

            if d_j <= TOL  # j is inside or on boundary
                push!(new_poly, p[j, :])
                if d_k > TOL  # k is outside -> add intersection
                    alpha = abs(d_j) / (abs(d_j) + abs(d_k))
                    intersection = p[j, :] .+ alpha .* (p[k, :] .- p[j, :])
                    push!(new_poly, intersection)
                end
            elseif d_k <= TOL  # j is outside, k is inside -> add intersection
                alpha = abs(d_j) / (abs(d_j) + abs(d_k))
                intersection = p[j, :] .+ alpha .* (p[k, :] .- p[j, :])
                push!(new_poly, intersection)
            end
        end

        if isempty(new_poly)
            return zeros(0, 2)
        end

        p = reduce(vcat, [v' for v in new_poly])
    end

    return p
end
