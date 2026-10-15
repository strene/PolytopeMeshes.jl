# Clipping of labeled convex polyhedra by half-spaces (replaces UPR's
# mirroredPebi3D/voronoi2mrstGrid3D, which relied on qhull).

"""
    ClipPolyhedron

Convex polyhedron with vertices `verts` and faces given as vertex loops
ordered counter-clockwise when seen from outside. `labels[f]` is the label of
face `f`: positive labels refer to neighboring sites, negative labels to
domain facets.
"""
struct ClipPolyhedron
    verts::Vector{Vec3}
    faces::Vector{Vector{Int}}
    labels::Vector{Int}
end

Base.isempty(P::ClipPolyhedron) = length(P.faces) < 4

"""
    clip(P, n, d, label, tol)

Clip the convex polyhedron `P` by the half-space `dot(n, x) <= d` (with `n` a
unit vector). The new face in the clipping plane gets `label`. Returns
`(Q, changed)`.
"""
function clip(P::ClipPolyhedron, n::Vec3, d::Float64, label::Int, tol::Float64)
    nv = length(P.verts)
    s = [dot(n, v) - d for v in P.verts]
    any(x -> x > tol, s) || return P, false
    all(x -> x > tol, s) && return ClipPolyhedron(Vec3[], Vector{Int}[], Int[]), true
    verts = copy(P.verts)
    edge_vertex = Dict{Tuple{Int, Int}, Int}()
    cap = Int[]
    faces = Vector{Int}[]
    labels = Int[]
    for (f, lab) in zip(P.faces, P.labels)
        m = length(f)
        loop = Int[]
        for k in 1:m
            a = f[k]
            b = f[k == m ? 1 : k + 1]
            sa, sb = s[a], s[b]
            if sa <= tol
                push!(loop, a)
                if sa >= -tol
                    push!(cap, a)
                end
            end
            if (sa > tol && sb < -tol) || (sa < -tol && sb > tol)
                key = a < b ? (a, b) : (b, a)
                id = get(edge_vertex, key, 0)
                if id == 0
                    t = sa / (sa - sb)
                    push!(verts, P.verts[a] + t * (P.verts[b] - P.verts[a]))
                    id = length(verts)
                    edge_vertex[key] = id
                end
                push!(loop, id)
                push!(cap, id)
            end
        end
        if length(loop) >= 3
            push!(faces, loop)
            push!(labels, lab)
        end
    end
    unique!(cap)
    if length(cap) >= 3
        # Order the cap vertices counter-clockwise around n
        c = sum(verts[i] for i in cap) / length(cap)
        u = verts[cap[1]] - c
        u = u - dot(u, n) * n
        if norm(u) < 1e-300
            u = abs(n[1]) < 0.9 ? cross(n, Vec3(1, 0, 0)) : cross(n, Vec3(0, 1, 0))
        end
        u = u / norm(u)
        w = cross(n, u)
        ang = [atan(dot(verts[i] - c, w), dot(verts[i] - c, u)) for i in cap]
        push!(faces, cap[sortperm(ang)])
        push!(labels, label)
    end
    # Remove unused vertices
    used = falses(length(verts))
    for f in faces, i in f
        used[i] = true
    end
    newid = cumsum(used)
    faces = [[newid[i] for i in f] for f in faces]
    return ClipPolyhedron(verts[used], faces, labels), true
end

function max_radius(P::ClipPolyhedron, x::Vec3)
    r = 0.0
    for v in P.verts
        r = max(r, sum(abs2, v - x))
    end
    return sqrt(r)
end

"""
    polytope_from_planes(normals, offsets)

Build the convex polyhedron `{x : dot(normals[k], x) <= offsets[k]}` by
clipping a large box. Faces are labeled `-k`.
"""
function polytope_from_planes(normals::Vector{Vec3}, offsets::Vector{Float64})
    # Large enough box: use the plane offsets as a scale
    R = 10 * (maximum(abs, offsets) + 1)
    box = ClipPolyhedron(
        [Vec3(-R, -R, -R), Vec3(R, -R, -R), Vec3(R, R, -R), Vec3(-R, R, -R),
         Vec3(-R, -R, R), Vec3(R, -R, R), Vec3(R, R, R), Vec3(-R, R, R)],
        [[1, 4, 3, 2], [5, 6, 7, 8], [1, 2, 6, 5], [3, 4, 8, 7], [1, 5, 8, 4], [2, 3, 7, 6]],
        fill(0, 6)
    )
    P = box
    tol = 1e-12 * R
    for (k, (n, d)) in enumerate(zip(normals, offsets))
        P, _ = clip(P, n, d, -k, tol)
    end
    any(==(0), P.labels) && throw(ArgumentError("The planes do not define a bounded domain"))
    return P
end
