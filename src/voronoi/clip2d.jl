# Clipping of labeled polygons by half-planes. A Voronoi cell restricted to a
# polygonal domain is obtained by clipping the domain by the bisectors of the
# site and its neighbors (UPR: clippedPebi2D / clipPolygon).

"""
    ClipPolygon

Polygon with counter-clockwise vertices `verts` where `labels[k]` labels the
edge from `verts[k]` to `verts[k+1]`. Positive labels refer to neighboring
sites, negative labels to domain boundary edges.
"""
struct ClipPolygon
    verts::Vector{Vec2}
    labels::Vector{Int}
end

Base.isempty(P::ClipPolygon) = length(P.verts) < 3

"""
    clip(P, n, d, label, tol)

Clip `P` by the half-plane `dot(n, x) <= d` (with `n` a unit vector). The new
edge along the clipping line gets `label`. Points within `tol` of the line are
considered to be on it. Returns `(Q, changed)`.
"""
function clip(P::ClipPolygon, n::Vec2, d::Float64, label::Int, tol::Float64)
    m = length(P.verts)
    s = [dot(n, v) - d for v in P.verts]
    any(x -> x > tol, s) || return P, false
    all(x -> x > tol, s) && return ClipPolygon(Vec2[], Int[]), true
    verts = Vec2[]
    labels = Int[]
    for k in 1:m
        kn = k == m ? 1 : k + 1
        a, b = s[k], s[kn]
        if a <= tol
            if b > tol
                if a < -tol
                    push!(verts, P.verts[k]); push!(labels, P.labels[k])
                    t = a / (a - b)
                    push!(verts, P.verts[k] + t * (P.verts[kn] - P.verts[k])); push!(labels, label)
                else
                    push!(verts, P.verts[k]); push!(labels, label)
                end
            else
                push!(verts, P.verts[k]); push!(labels, P.labels[k])
            end
        elseif b < -tol
            t = a / (a - b)
            push!(verts, P.verts[k] + t * (P.verts[kn] - P.verts[k])); push!(labels, P.labels[k])
        end
    end
    return ClipPolygon(verts, labels), true
end

function max_radius(P::ClipPolygon, x::Vec2)
    r = 0.0
    for v in P.verts
        r = max(r, sum(abs2, v - x))
    end
    return sqrt(r)
end
