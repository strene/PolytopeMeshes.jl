# Basic 2D polygon and segment geometry.

cross2(a, b) = a[1] * b[2] - a[2] * b[1]

"""
    signed_area(poly)

Signed area of a closed polygon given by its vertices (positive if the
vertices are ordered counter-clockwise).
"""
function signed_area(poly::AbstractVector{<:SVector{2}})
    A = 0.0
    n = length(poly)
    for i in 1:n
        A += cross2(poly[i], poly[mod1(i + 1, n)])
    end
    return A / 2
end

"""
    inpolygon(x, poly)

Test if the point `x` lies inside (or on the boundary of) the polygon `poly`.
"""
function inpolygon(x::SVector{2}, poly::AbstractVector{<:SVector{2}}; tol = 0.0)
    n = length(poly)
    inside = false
    for i in 1:n
        a = poly[i]
        b = poly[mod1(i + 1, n)]
        if tol > 0 && point_segment_distance(x, a, b) <= tol
            return true
        end
        if (a[2] > x[2]) != (b[2] > x[2])
            xc = a[1] + (x[2] - a[2]) * (b[1] - a[1]) / (b[2] - a[2])
            if x[1] < xc
                inside = !inside
            end
        end
    end
    return inside
end

function point_segment_distance(x, a, b)
    v = b - a
    L2 = dot(v, v)
    t = L2 > 0 ? clamp(dot(x - a, v) / L2, 0.0, 1.0) : 0.0
    return norm(x - (a + t * v))
end

"""
    signed_distance_polygon(x, poly)

Signed distance from `x` to the boundary of `poly`, negative inside (the
`dpoly` function of DistMesh).
"""
function signed_distance_polygon(x::SVector{2}, poly::AbstractVector{<:SVector{2}})
    n = length(poly)
    d = Inf
    for i in 1:n
        d = min(d, point_segment_distance(x, poly[i], poly[mod1(i + 1, n)]))
    end
    return inpolygon(x, poly) ? -d : d
end

"""
    signed_distance_rectangle(x, lo, hi)

Signed distance (approximate at the corners, as in DistMesh's `drectangle`)
from `x` to the rectangle `[lo, hi]`.
"""
function signed_distance_rectangle(x::SVector{2}, lo, hi)
    return -min(min(min(-lo[2] + x[2], hi[2] - x[2]), -lo[1] + x[1]), hi[1] - x[1])
end

"""
    segment_intersection(p1, p2, q1, q2; tol)

Intersection of the segments `p1-p2` and `q1-q2`. Returns `(t, u)` such that
the intersection is `p1 + t (p2 - p1) = q1 + u (q2 - q1)`, or `nothing` if
the segments do not intersect (parallel segments never intersect).
"""
function segment_intersection(p1, p2, q1, q2; tol = 1e-12)
    r = p2 - p1
    s = q2 - q1
    den = cross2(r, s)
    scale = norm(r) * norm(s)
    if abs(den) <= 1e-14 * scale
        return nothing
    end
    qp = q1 - p1
    t = cross2(qp, s) / den
    u = cross2(qp, r) / den
    if -tol <= t <= 1 + tol && -tol <= u <= 1 + tol
        return (clamp(t, 0.0, 1.0), clamp(u, 0.0, 1.0))
    end
    return nothing
end

"""
    clip_segment_convex(a, b, poly)

Clip the segment `a-b` against the convex counter-clockwise polygon `poly`.
Returns the clipped end points, or `nothing` if the segment is outside.
"""
function clip_segment_convex(a::SVector{2}, b::SVector{2}, poly::AbstractVector{<:SVector{2}}; tol = 1e-12)
    t0, t1 = 0.0, 1.0
    d = b - a
    n = length(poly)
    for i in 1:n
        p = poly[i]
        q = poly[mod1(i + 1, n)]
        e = q - p
        # Outward normal of a counter-clockwise polygon edge
        nrm = SVector(e[2], -e[1])
        num = dot(nrm, p - a)
        den = dot(nrm, d)
        if abs(den) < 1e-14 * norm(nrm) * max(norm(d), eps())
            if num < -tol * norm(nrm)
                return nothing
            end
        else
            t = num / den
            if den < 0
                t0 = max(t0, t)
            else
                t1 = min(t1, t)
            end
            if t0 > t1
                return nothing
            end
        end
    end
    return a + t0 * d, a + t1 * d
end
