# 2D grid generators (UPR: clippedPebi2D, pebiGrid2D, compositePebiGrid2D).

"""
    clipped_voronoi_grid_2d(sites, boundary)

Voronoi grid of `sites` clipped against the polygon `boundary` (UPR:
`clippedPebi2D`).
"""
clipped_voronoi_grid_2d(sites, boundary; kwarg...) = voronoi_mesh(sites, PolygonDomain(boundary); kwarg...)

function _domain_2d(pdims, polygon)
    if isnothing(polygon)
        isnothing(pdims) && throw(ArgumentError("Either pdims or polygon must be given"))
        length(pdims) == 2 || throw(ArgumentError("pdims must have 2 elements"))
        all(pdims .> 0) || throw(ArgumentError("pdims must be positive"))
        return box_domain((0.0, 0.0), (pdims[1], pdims[2]))
    end
    return PolygonDomain(polygon)
end

# Sites along cell constraints for the 2D generators. Lines are split at their
# intersections and at intersections with face constraints; point constraints
# are passed through.
function _cell_constraint_sites_2d(cell_constraints, face_constraints; cc_grid_size, interpolate_cc, cc_rho,
        circle_factor, protection_layer, protection_distance, fault_length)
    cc_all = to_polylines(cell_constraints, 2)
    fc_all = to_polylines(face_constraints, 2)
    ncc = length(cc_all)
    interpolate_cc = _per_line(interpolate_cc, ncc)
    protection_distance = _per_line(protection_distance, ncc)
    line_ids = findall(l -> length(l) > 1, cc_all)
    point_ids = findall(l -> length(l) == 1, cc_all)
    split = split_at_intersections(cc_all[line_ids], fc_all)
    origin = line_ids[split.origin]
    # Stop the cell constraint sites before crossings with face constraints so
    # that the face constraint sites can represent the crossing.
    bisect = fault_length / 2
    fault_offset = sqrt((circle_factor * fault_length)^2 - bisect^2)
    fac = 1.0 + fault_offset / cc_grid_size
    se = [(c in (2, 3) ? fac : 0.0, c in (1, 3) ? fac : 0.0) for c in split.other_cut]
    lines = vcat(split.lines, cc_all[point_ids])
    origin = vcat(origin, point_ids)
    se = vcat(se, fill((0.0, 0.0), length(point_ids)))
    cut = vcat(split.self_cut, zeros(Int, length(point_ids)))
    return line_sites_2d(lines, cc_grid_size;
        rho = cc_rho,
        interpolate = interpolate_cc[origin],
        se = se,
        cut = cut,
        protection_layer = protection_layer,
        protection_distance = protection_distance[origin],
        line_index = origin
    )
end

function _voronoi_mesh_2d(F, W, other_sites, domain; kwarg...)
    nF = length(F.sites)
    nW = length(W.sites)
    sites = vcat(F.sites, W.sites, other_sites)
    tags = zeros(Int, length(sites))
    tags[1:nF] .= crossing_cell_tags(F)
    tags[(nF + 1):(nF + nW)] .= W.line_index
    face_tag = (i, j) -> (i <= nF && j <= nF) ? surface_face_tag(F, i, j) : 0
    return voronoi_mesh(sites, domain; cell_tags = tags, face_tag = face_tag, kwarg...)
end

"""
    voronoi_grid_2d(grid_size, pdims; kwarg...)

Fully unstructured 2D Voronoi grid that conforms to face constraints (faults,
fractures) and cell constraints (wells, points) (UPR: `pebiGrid2D`). The
background sites are distributed with DistMesh.

The domain is the rectangle `[0, pdims[1]] × [0, pdims[2]]` unless `polygon`
is given.

# Keyword arguments
- `polygon = nothing`: vertices of a polygonal domain (replaces `pdims`).
- `cell_constraints = []`: polylines; the cells are placed along them. A
  single point gives a point constraint.
- `cc_factor = 1`: relative size of the constrained cells.
- `interpolate_cc = false`: per constraint; if `false` each polyline segment
  is represented by at least one cell.
- `cc_refinement = false`: refine the background grid around cell
  constraints.
- `cc_eps = 0.25*maximum(pdims)`: refinement transition length.
- `cc_rho = x -> 1`: relative spacing along the cell constraints.
- `protection_layer = false`: add sites on both sides of the cell
  constraints.
- `protection_distance = grid_size/10`: distance to the protection sites.
- `face_constraints = []`: polylines the faces should conform to.
- `fc_factor = 1`: relative size of the face constraint cells.
- `interpolate_fc = false`: per constraint, as `interpolate_cc`.
- `circle_factor = 0.6`: see [`surface_sites_2d`](@ref).
- `fc_rho = x -> 1`: relative spacing along the face constraints.
- `fc_refinement = false`, `fc_eps`: refinement around face constraints.
- `sufficient_condition = true`: remove background sites inside the face
  constraint circles (guarantees conformity). If `false`, only sites closer
  than the face constraint spacing are removed.
- `triangulation = DistMeshBackend()`: backend used to distribute the
  background sites, see [`TriangulationBackend`](@ref).
- `rng`: random number generator passed to the backend (seeded by default
  for reproducible grids).

Returns a [`PolytopalMesh`](@ref). Cells generated from cell constraints are
tagged with the constraint index in `cell_tags`, faces on face constraints
with the constraint index in `face_tags`.
"""
function voronoi_grid_2d(grid_size, pdims = nothing;
        polygon = nothing,
        cell_constraints = Vector{Vec2}[],
        cc_factor = 1.0,
        interpolate_cc = false,
        cc_refinement = false,
        cc_rho = x -> 1.0,
        cc_eps = nothing,
        protection_layer = false,
        protection_distance = grid_size / 10,
        face_constraints = Vector{Vec2}[],
        fc_factor = 1.0,
        interpolate_fc = false,
        circle_factor = 0.6,
        fc_refinement = false,
        fc_rho = x -> 1.0,
        fc_eps = nothing,
        sufficient_condition = true,
        triangulation::TriangulationBackend = default_triangulation_backend(),
        rng = MersenneTwister(1234),
        kwarg...
    )
    grid_size > 0 || throw(ArgumentError("grid_size must be positive"))
    0.5 < circle_factor < 1 || throw(ArgumentError("circle_factor must be in (0.5, 1)"))
    domain = _domain_2d(pdims, polygon)
    lo, hi = domain_bounding_box(domain)
    L = maximum(hi - lo)
    cc_eps = something(cc_eps, 0.25 * L)
    fc_eps = something(fc_eps, 0.25 * L)
    cc_size = grid_size * cc_factor
    fc_size = grid_size * fc_factor
    cc_rho_abs = x -> cc_size * _eval_size(cc_rho, x)
    fault_length = cc_refinement ? cc_size * fc_factor * 1.2 : fc_size
    W = _cell_constraint_sites_2d(cell_constraints, face_constraints; cc_grid_size = cc_size,
        interpolate_cc = interpolate_cc, cc_rho = cc_rho_abs, circle_factor = circle_factor,
        protection_layer = protection_layer, protection_distance = protection_distance,
        fault_length = fault_length)
    # Background size near cell constraints (absolute)
    h_cc = x -> Float64(grid_size)
    if cc_refinement && !isempty(W.sites)
        wtree = KDTree(W.sites)
        h_cc = x -> min(grid_size, 1.2 * cc_size * exp(nn(wtree, x)[2] / cc_eps))
    end
    if cc_refinement
        fault_spacing = x -> fc_factor * h_cc(x) * _eval_size(fc_rho, x)
    else
        fault_spacing = x -> fc_size * _eval_size(fc_rho, x)
    end
    F = surface_sites_2d(face_constraints, fc_size; fh = fault_spacing, circle_factor = circle_factor,
        interpolate = interpolate_fc, cell_constraints = cell_constraints)
    tips = filter(x -> inside(domain, x), F.tips)
    # Background size function (absolute)
    h_fc = x -> Float64(grid_size)
    if fc_refinement && !isempty(F.sites)
        ftree = KDTree(F.sites)
        h_fc = x -> min(grid_size, 1.2 * fc_size * exp(nn(ftree, x)[2] / fc_eps))
    end
    h0 = min(cc_refinement ? cc_size : grid_size, fc_refinement ? fc_size : grid_size)
    h0 = min(h0, grid_size)
    hfun = x -> min(h_cc(x), h_fc(x))
    fixed = vcat(F.sites, W.sites, W.protection_sites, tips, domain.vertices)
    nconstr = length(F.sites) + length(W.sites)
    pts, _, _ = _triangulate(triangulation, domain.vertices, hfun, fixed; rng = rng, min_size = h0)
    background = pts[(nconstr + 1):end]
    if sufficient_condition
        background, _ = surface_sufficient_condition(background, F.centers, F.radii)
    else
        background, _ = remove_conflict_points(background, F.sites, F.spacing)
    end
    return _voronoi_mesh_2d(F, W, background, domain; kwarg...)
end

"""
    composite_voronoi_grid_2d(cell_size, pdims; kwarg...)

2D Voronoi grid with a Cartesian background grid that is refined around cell
constraints and conforms to face and cell constraints (UPR:
`compositePebiGrid2D`).

`cell_size = (dx, dy)` is the size of the background cells. The domain is the
rectangle `[0, pdims[1]] × [0, pdims[2]]` unless `polygon` is given.

# Keyword arguments
- `polygon`, `cell_constraints`, `cc_factor`, `interpolate_cc`, `cc_rho`,
  `protection_layer`, `protection_distance`, `face_constraints`,
  `fc_factor`, `interpolate_fc`, `circle_factor`: see [`voronoi_grid_2d`](@ref).
- `mlqt_max_level = 0`: number of quadtree refinement levels around the cell
  constraints.
- `mlqt_level_steps = -1`: refinement distance per level, see
  [`mlqt`](@ref).
"""
function composite_voronoi_grid_2d(cell_size, pdims = nothing;
        polygon = nothing,
        cell_constraints = Vector{Vec2}[],
        cc_factor = 1.0,
        interpolate_cc = false,
        mlqt_max_level = 0,
        mlqt_level_steps = -1,
        cc_rho = x -> 1.0,
        face_constraints = Vector{Vec2}[],
        fc_factor = 1.0,
        interpolate_fc = false,
        circle_factor = 0.6,
        protection_layer = false,
        protection_distance = norm(cell_size) / 10,
        kwarg...
    )
    length(cell_size) == 2 || throw(ArgumentError("cell_size must have 2 elements"))
    all(cell_size .> 0) || throw(ArgumentError("cell_size must be positive"))
    mlqt_max_level >= 0 || throw(ArgumentError("mlqt_max_level must be non-negative"))
    0.5 < circle_factor < 1 || throw(ArgumentError("circle_factor must be in (0.5, 1)"))
    domain = _domain_2d(pdims, polygon)
    cc_size = minimum(cell_size) * cc_factor
    fc_size = minimum(cell_size) * fc_factor
    cc_rho_abs = x -> cc_size * _eval_size(cc_rho, x)
    W = _cell_constraint_sites_2d(cell_constraints, face_constraints; cc_grid_size = cc_size,
        interpolate_cc = interpolate_cc, cc_rho = cc_rho_abs, circle_factor = circle_factor,
        protection_layer = protection_layer, protection_distance = protection_distance,
        fault_length = fc_size)
    F = surface_sites_2d(face_constraints, fc_size; circle_factor = circle_factor,
        interpolate = interpolate_fc, cell_constraints = cell_constraints)
    # Cartesian background sites (including the domain boundary)
    lo, hi = domain_bounding_box(domain)
    nx = ceil(Int, (hi[1] - lo[1]) / cell_size[1] - 1e-10)
    ny = ceil(Int, (hi[2] - lo[2]) / cell_size[2] - 1e-10)
    dx = (hi[1] - lo[1]) / nx
    dy = (hi[2] - lo[2]) / ny
    background = Vec2[]
    for j in 0:ny, i in 0:nx
        x = Vec2(lo[1] + i * dx, lo[2] + j * dy)
        if isnothing(polygon) || inside(domain, x)
            push!(background, x)
        end
    end
    if !isempty(W.sites) && mlqt_max_level > 0
        tree = KDTree(W.sites)
        refined = Vec2[]
        for x in background
            for (c, _) in mlqt(x, W.sites, (dx, dy); max_level = mlqt_max_level, dist_tol = mlqt_level_steps, tree = tree)
                push!(refined, c)
            end
        end
        background = refined
    end
    tips = filter(x -> inside(domain, x), F.tips)
    background, _ = remove_conflict_points(background, W.sites, W.spacing)
    background, _ = remove_conflict_points(background, W.protection_sites, W.protection_spacing)
    background, _ = remove_conflict_points(background, F.sites, F.spacing)
    background, _ = remove_conflict_points(background, F.centers, F.radii)
    other = vcat(W.protection_sites, tips, background)
    return _voronoi_mesh_2d(F, W, other, domain; kwarg...)
end
