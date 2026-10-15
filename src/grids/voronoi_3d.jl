# 3D grid generators (UPR: mirroredPebi3D, compositePebiGrid3D).

"""
    mirrored_voronoi_grid_3d(sites, boundary_points)

Voronoi grid of `sites` restricted to the convex hull of
`boundary_points` (UPR: `mirroredPebi3D`).
"""
mirrored_voronoi_grid_3d(sites, boundary_points; kwarg...) = voronoi_mesh(sites, convex_hull_domain(boundary_points); kwarg...)

"""
    composite_voronoi_grid_3d(celldim, pdims; kwarg...)

3D Voronoi grid with a Cartesian background grid of `celldim = (nx, ny, nz)`
cells in the box `[0, pdims]`, refined around cell constraints, that conforms
to triangulated face constraints and to cell constraints (UPR:
`compositePebiGrid3D`).

# Keyword arguments
- `cell_constraints = []`: 3D polylines along which cells are placed.
- `cc_rho`: spacing of the sites along the cell constraints (number, function
  of position, or one per constraint). Defaults to the smallest background
  cell size.
- `mlqt_max_level = 0`, `mlqt_level_steps = -1`: octree refinement around
  the cell constraints, see [`mlqt`](@ref).
- `face_constraints = []`: triangulated surfaces `(points, triangles)`, see
  [`surface_sites_3d`](@ref) and [`surface_triangulation`](@ref).
- `fc_rho`: ball radius at the surface vertices (number, function, or one per
  surface). Defaults to the smallest background cell size.

Tip sites beyond free surface edges inside the domain make the conforming
faces stop at the edges.

Cells generated from cell constraints are tagged with the constraint index in
`cell_tags`, faces on the surfaces with the surface index in `face_tags`.
"""
function composite_voronoi_grid_3d(celldim, pdims;
        cell_constraints = Vector{Vec3}[],
        cc_rho = nothing,
        mlqt_max_level = 0,
        mlqt_level_steps = -1,
        face_constraints = [],
        fc_rho = nothing,
        kwarg...
    )
    length(celldim) == 3 || throw(ArgumentError("celldim must have 3 elements"))
    length(pdims) == 3 || throw(ArgumentError("pdims must have 3 elements"))
    all(celldim .> 0) || throw(ArgumentError("celldim must be positive"))
    all(pdims .> 0) || throw(ArgumentError("pdims must be positive"))
    mlqt_max_level >= 0 || throw(ArgumentError("mlqt_max_level must be non-negative"))
    h = Vec3(pdims ./ celldim)
    cc_rho = something(cc_rho, minimum(h))
    fc_rho = something(fc_rho, minimum(h))
    domain = box_domain((0.0, 0.0, 0.0), Tuple(Float64.(pdims)))
    F = surface_sites_3d(face_constraints, fc_rho; domain = domain)
    F = remove_surface_conflict_sites_3d(F)
    W = line_sites_3d(cell_constraints, cc_rho)
    # Cartesian background sites at the cell centers
    background = Vec3[]
    for k in 1:celldim[3], j in 1:celldim[2], i in 1:celldim[1]
        push!(background, Vec3((i - 0.5) * h[1], (j - 0.5) * h[2], (k - 0.5) * h[3]))
    end
    if !isempty(W.sites) && mlqt_max_level > 0
        tree = KDTree(W.sites)
        refined = Vec3[]
        for x in background
            for (c, _) in mlqt(x, W.sites, h; max_level = mlqt_max_level, dist_tol = mlqt_level_steps, tree = tree)
                push!(refined, c)
            end
        end
        background = refined
    end
    background, _ = line_sufficient_condition(background, W.paths)
    background, _ = surface_sufficient_condition(background, F.centers, F.radii)
    wsites, wrem = surface_sufficient_condition(W.sites, F.centers, F.radii)
    windex = W.line_index[.!wrem]
    nF = length(F.sites)
    sites = vcat(F.sites, wsites, F.tips, background)
    tags = zeros(Int, length(sites))
    tags[(nF + 1):(nF + length(wsites))] .= windex
    face_tag = (i, j) -> (i <= nF && j <= nF) ? surface_face_tag(F, i, j) : 0
    return voronoi_mesh(sites, domain; cell_tags = tags, face_tag = face_tag, kwarg...)
end

"""
    surface_triangulation(polygon, h; kwarg...)

Triangulate the planar polygon `polygon` (vertices in 3D, or in 2D) with
approximately equilateral triangles of size `h`. Returns `(points,
triangles)` for use with [`surface_sites_3d`](@ref) or
[`composite_voronoi_grid_3d`](@ref).

The triangulation is created by `triangulation` (default
[`DistMeshBackend`](@ref)), see [`TriangulationBackend`](@ref).
"""
function surface_triangulation(polygon, h; triangulation::TriangulationBackend = default_triangulation_backend(),
        rng = MersenneTwister(1234))
    D = polygon isa AbstractMatrix ? size(polygon, 2) : length(first(polygon))
    if D == 2
        poly2 = to_points(polygon, 2)
        R, c = SMatrix{3, 3, Float64}(I), zero(Vec3)
    else
        poly3 = to_points(polygon, 3)
        R, c = plane_frame(poly3)
        poly2 = [to_plane(R, c, x) for x in poly3]
    end
    if signed_area(poly2) < 0
        poly2 = reverse(poly2)
    end
    # Fixed points along the polygon edges
    fixed = Vec2[]
    n = length(poly2)
    for i in 1:n
        a, b = poly2[i], poly2[mod1(i + 1, n)]
        m = max(ceil(Int, norm(b - a) / h - 1e-10), 1)
        for k in 0:(m - 1)
            push!(fixed, a + k / m * (b - a))
        end
    end
    pts, triangles, _ = _triangulate(triangulation, poly2, x -> Float64(h), fixed; rng = rng, min_size = h)
    isempty(triangles) && error("Backend $(typeof(triangulation)) returned no triangles")
    # Remove points not used by any triangle
    used = falses(length(pts))
    for t in triangles, i in t
        used[i] = true
    end
    newid = cumsum(used)
    triangles = [(newid[i], newid[j], newid[k]) for (i, j, k) in triangles]
    pts = pts[used]
    points = D == 2 ? [Vec3(p[1], p[2], 0.0) for p in pts] : [from_plane(R, c, p) for p in pts]
    return (points = points, triangles = triangles)
end
