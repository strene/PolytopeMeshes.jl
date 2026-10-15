# Hierarchical 3D gridding of planar fractures where the 3D grid conforms
# exactly to 2D grids of the fractures, which in turn conform to 1D grids of
# the fracture intersections (UPR: ellipticSurface3D, surfaceIntersections3D,
# lineGrid3D, surfaceGrid3D, volumeGrid3D; Berge et al., Comput Geosci 2019).

"""
    elliptic_surface_3d(center, major_axis, minor_axis, major_axis_angle, strike_angle, dip_angle; num_points = 16)

Polygon approximating an elliptic surface (UPR: `ellipticSurface3D`, adapted
from PorePy). The ellipse with half axes `major_axis` and `minor_axis` is
rotated by `major_axis_angle` around the z-axis, then by `dip_angle` around
the strike direction, which is given by `strike_angle` measured from the
x-axis. Angles are in radians.
"""
function elliptic_surface_3d(center, major_axis, minor_axis, major_axis_angle, strike_angle, dip_angle; num_points = 16)
    angs = range(0, 2pi, length = num_points + 1)[1:(end - 1)]
    ref = [Vec3(major_axis * cos(a), minor_axis * sin(a), 0.0) for a in angs]
    Rmajor = rotation_matrix(major_axis_angle, Vec3(0, 0, 1))
    strike_dir = rotation_matrix(strike_angle, Vec3(0, 0, 1)) * Vec3(1, 0, 0)
    Rdip = rotation_matrix(dip_angle, strike_dir)
    c = Vec3(center)
    return [c + Rdip * (Rmajor * p) for p in ref]
end

"""
    polygon_intersection_3d(poly1, poly2; tol = 1e-6)

Intersection segment of two convex planar polygons in 3D, or `nothing`.
Parallel polygons are not treated.
"""
function polygon_intersection_3d(poly1::AbstractVector{Vec3}, poly2::AbstractVector{Vec3}; tol = 1e-6)
    R, c = plane_frame(poly1)
    p1 = [to_plane(R, c, x) for x in poly1]
    if signed_area(p1) < 0
        p1 = reverse(p1)
    end
    q = [R * (x - c) for x in poly2]
    z = [x[3] for x in q]
    if all(z .> tol) || all(z .< -tol) || all(abs.(z) .< tol)
        return nothing
    end
    AB = Vec2[]
    n = length(q)
    for e in 1:n
        v1, v2 = q[e], q[mod1(e + 1, n)]
        s1, s2 = sign(v1[3]), sign(v2[3])
        s1 == s2 && continue
        d1, d2 = abs(v1[3]), abs(v2[3])
        t = d1 / (d1 + d2)
        I = v1 + (v2 - v1) * t
        push!(AB, Vec2(I[1], I[2]))
    end
    AB, _ = unique_points(AB; tol = 1e-10)
    length(AB) == 2 || return nothing
    seg = clip_segment_convex(AB[1], AB[2], p1)
    isnothing(seg) && return nothing
    a, b = seg
    norm(b - a) < tol && return nothing
    return from_plane(R, c, a), from_plane(R, c, b)
end

"""
    surface_intersections_3d(surfaces)

Intersections between convex planar polygons (UPR: `surfaceIntersections3D`).
Returns a vector of named tuples `(points = (a, b), surfaces = (f1, f2))`.
Intersection lines that cross each other in a surface are split at the
crossing.
"""
function surface_intersections_3d(surfaces)
    surfaces = [to_points(s, 3) for s in surfaces]
    ints = NamedTuple{(:points, :surfaces), Tuple{Tuple{Vec3, Vec3}, Tuple{Int, Int}}}[]
    for f1 in eachindex(surfaces), f2 in (f1 + 1):length(surfaces)
        seg = polygon_intersection_3d(surfaces[f1], surfaces[f2])
        isnothing(seg) && continue
        push!(ints, (points = seg, surfaces = (f1, f2)))
    end
    # Split intersection lines at their crossings within each surface
    for f in eachindex(surfaces)
        idx = findall(x -> f in x.surfaces, ints)
        length(idx) < 2 && continue
        R, c = plane_frame(surfaces[f])
        lines = [[to_plane(R, c, ints[j].points[1]), to_plane(R, c, ints[j].points[2])] for j in idx]
        split = split_at_intersections(lines)
        new = [(points = (from_plane(R, c, l[1]), from_plane(R, c, l[end])), surfaces = ints[idx[o]].surfaces)
            for (l, o) in zip(split.lines, split.origin)]
        ints = vcat(ints[setdiff(eachindex(ints), idx)], new)
    end
    return ints
end

"""
    line_grid_3d(intersections, ds)

1D grids of the intersection lines (UPR: `lineGrid3D`). For each intersection
returns a named tuple `(sites, vertices)` with the 1D cell centers and the
`n+1` cell vertices. Intersections shorter than `1.5*ds` get no cells.
"""
function line_grid_3d(intersections, ds)
    out = NamedTuple{(:sites, :vertices), Tuple{Vector{Vec3}, Vector{Vec3}}}[]
    for int in intersections
        p1, p2 = int.points
        v = p2 - p1
        L = norm(v)
        v = v / L
        if L < 1.5 * ds
            push!(out, (sites = Vec3[], vertices = Vec3[]))
            continue
        end
        n = round(Int, L / ds)
        t = collect(range(ds / 2, L - ds / 2, length = n))
        tv = vcat(0.0, (t[1:(end - 1)] .+ t[2:end]) ./ 2, L)
        push!(out, (sites = [p1 + v * s for s in t], vertices = [p1 + v * s for s in tv]))
    end
    return out
end

"""
    SurfaceGrid

2D grid of a planar surface embedded in 3D, as returned by
[`surface_grid_3d`](@ref).

# Fields
- `mesh`: 2D [`PolytopalMesh`](@ref) in the local coordinates of the plane.
  Faces conforming to intersection lines are tagged with the intersection
  index.
- `rotation`, `center`: map local coordinates `y` to 3D by
  `rotation' * [y; 0] + center`.
- `sites`: cell sites in 3D.
- `nodes`: node coordinates in 3D.
- `normal`: unit normal of the plane.
"""
struct SurfaceGrid
    mesh::PolytopalMesh{2}
    rotation::SMatrix{3, 3, Float64, 9}
    center::Vec3
    sites::Vector{Vec3}
    nodes::Vector{Vec3}
    normal::Vec3
end

Base.show(io::IO, g::SurfaceGrid) = print(io, "SurfaceGrid with $(number_of_cells(g.mesh)) cells")

function _suf_cond_from_1d(pts, grids_1, gamma)
    for (c, v) in grids_1
        cc = vcat(c, [c[end]])
        R = [sqrt(sum(abs2, cc[k] - v[k]) + gamma^2) for k in eachindex(v)]
        pts, _ = remove_conflict_points(pts, v, R)
    end
    return pts
end

"""
    surface_grid_3d(surfaces, grids_1, intersections, ds, gamma; triangulation, rng)

2D Voronoi grids of planar polygonal surfaces conforming to the 1D grids of
their intersections (UPR: `surfaceGrid3D`). `ds` is the 2D cell size and
`gamma` the distance from the intersection lines to the 2D sites next to
them (reduce it if there are problems at intersections). The background
sites are distributed by `triangulation` (default [`DistMeshBackend`](@ref)).
Returns a vector of [`SurfaceGrid`](@ref).
"""
function surface_grid_3d(surfaces, grids_1, intersections, ds, gamma;
        triangulation::TriangulationBackend = default_triangulation_backend(), rng = MersenneTwister(1234))
    surfaces = [to_points(s, 3) for s in surfaces]
    out = SurfaceGrid[]
    for f in eachindex(surfaces)
        R, c = plane_frame(surfaces[f])
        f_xy = [to_plane(R, c, x) for x in surfaces[f]]
        domain = PolygonDomain(f_xy)
        internal = Vec2[]
        internal_int = Int[]
        internal_id = Int[]
        tips = Vec2[]
        g1 = Tuple{Vector{Vec2}, Vector{Vec2}}[]
        for (j, int) in enumerate(intersections)
            f in int.surfaces || continue
            pts1 = grids_1[j].sites
            isempty(pts1) && continue
            p = [to_plane(R, c, x) for x in pts1]
            v = [to_plane(R, c, x) for x in grids_1[j].vertices]
            push!(g1, (p, v))
            t = p[end] - p[1]
            t = t / norm(t)
            nrm = Vec2(t[2], -t[1])
            for (k, x) in enumerate(p)
                push!(internal, x + nrm * gamma); push!(internal_int, j); push!(internal_id, k)
            end
            for (k, x) in enumerate(p)
                push!(internal, x - nrm * gamma); push!(internal_int, j); push!(internal_id, -k)
            end
            Rs = sqrt(sum(abs2, v[1] - p[1]) + gamma^2)
            Re = sqrt(sum(abs2, v[end] - p[end]) + gamma^2)
            push!(tips, v[1] - t * Rs)
            push!(tips, v[end] + t * Re)
        end
        tips = filter(x -> inside(domain, x), tips)
        tips = _suf_cond_from_1d(tips, g1, gamma * (1 - 1e-6))
        fixed = vcat(internal, tips, domain.vertices)
        pts, _, _ = _triangulate(triangulation, domain.vertices, x -> Float64(ds), fixed; rng = rng, min_size = ds)
        res = _suf_cond_from_1d(pts[(length(internal) + 1):end], g1, gamma * (1 - 1e-6))
        sites = vcat(internal, res)
        ni = length(internal)
        face_tag = (a, b) -> begin
            if a <= ni && b <= ni && internal_int[a] == internal_int[b] && internal_id[a] == -internal_id[b]
                return internal_int[a]
            end
            return 0
        end
        mesh = voronoi_mesh(sites, domain; face_tag = face_tag)
        sites3 = [from_plane(R, c, x) for x in mesh.cell_sites]
        nodes3 = [from_plane(R, c, x) for x in mesh.node_points]
        push!(out, SurfaceGrid(mesh, R, c, sites3, nodes3, Vec3(R[3, 1], R[3, 2], R[3, 3])))
    end
    return out
end

"""
    surface_sufficient_condition_from_grid(pts, grids_2, gamma)

Remove the points inside the balls centered at the nodes of the 2D surface
grids with radius `sqrt(kappa^2 + gamma^2)`, where `kappa` is the distance from
the node to the sites of the adjacent cells (UPR: `surfaceSufCondFromGrid3D`).
"""
function surface_sufficient_condition_from_grid(pts::AbstractVector{Vec3}, grids_2, gamma)
    removed = falses(length(pts))
    for g in grids_2
        m = g.mesh
        kappa2 = zeros(number_of_nodes(m))
        for f in 1:number_of_faces(m)
            for c in m.face_neighbors[f]
                c == 0 && continue
                for n in face_nodes(m, f)
                    kappa2[n] = max(kappa2[n], sum(abs2, g.nodes[n] - g.sites[c]))
                end
            end
        end
        R = sqrt.(kappa2 .+ gamma^2)
        _, rem = surface_sufficient_condition(pts, g.nodes, R)
        removed .|= rem
    end
    return pts[.!removed], removed
end

"""
    volume_grid_3d(pdims, surfaces, grids_2, ds, gamma)

3D Voronoi grid of the box `[0, pdims]` conforming exactly to the 2D surface
grids `grids_2` (UPR: `volumeGrid3D`). The sites next to the surfaces are
placed at distance `gamma` on both sides of the 2D sites, and the background
sites form a Cartesian lattice of size about `1.2*ds`. Faces on the surfaces
are tagged with the surface index.
"""
function volume_grid_3d(pdims, surfaces, grids_2, ds, gamma; kwarg...)
    internal = Vec3[]
    gen = Tuple{Int, Int, Int}[]
    for (f, g) in enumerate(grids_2)
        for (k, s) in enumerate(g.sites)
            push!(internal, s + g.normal * gamma); push!(gen, (f, k, 1))
        end
        for (k, s) in enumerate(g.sites)
            push!(internal, s - g.normal * gamma); push!(gen, (f, k, -1))
        end
    end
    h = 1.2 * ds
    n = [max(round(Int, pdims[d] / h), 1) for d in 1:3]
    background = Vec3[]
    for k in 1:n[3], j in 1:n[2], i in 1:n[1]
        push!(background, Vec3((i - 0.5) * pdims[1] / n[1], (j - 0.5) * pdims[2] / n[2], (k - 0.5) * pdims[3] / n[3]))
    end
    background, _ = surface_sufficient_condition_from_grid(background, grids_2, gamma)
    sites = vcat(internal, background)
    ni = length(internal)
    face_tag = (a, b) -> begin
        if a <= ni && b <= ni
            fa, ka, sa = gen[a]
            fb, kb, sb = gen[b]
            if fa == fb && ka == kb && sa == -sb
                return fa
            end
        end
        return 0
    end
    domain = box_domain((0.0, 0.0, 0.0), Tuple(Float64.(pdims)))
    return voronoi_mesh(sites, domain; face_tag = face_tag, kwarg...)
end
