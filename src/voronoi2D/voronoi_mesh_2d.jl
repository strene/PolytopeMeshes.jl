"""
Construct a 2D Voronoi mesh using DistMesh for background grid generation.
Ported from MRST's pebiGrid2D (UPR module), renamed following the package's
naming convention (pebi → voronoi).
"""

"""
    voronoi_mesh_2d(res_grid_size, pdims; kwargs...)

Construct a 2D Voronoi mesh with adaptive background grid generation using
DistMesh. Supports cell constraints (wells) and face constraints (faults)
with optional refinement.

This is the DistMesh-based alternative to [`composite_voronoi_mesh_2d`](@ref),
which uses a Cartesian background grid. The DistMesh approach produces more
uniform cell sizes and supports density-based refinement around constraints.

# Arguments
- `res_grid_size`: Size of the reservoir grid cells (scalar, in meters).
- `pdims`: `[xmax, ymax]` physical dimensions of the domain.

# Keyword Arguments
- `cell_constraints`: Vector of n×2 matrices defining well/cell constraint paths
  (default: empty)
- `cc_factor`: Relative grid size for cell constraints compared to reservoir
  cells (default: 1.0). If `cc_factor=0.5`, constrained cells are about half
  the reservoir cell size.
- `interpolate_cc`: Whether to interpolate cell constraints (default: false).
  Can be a single `Bool` or a `Vector{Bool}` with one entry per constraint.
- `cc_refinement`: Enable refinement around cell constraints (default: false)
- `cc_eps`: Refinement transition distance around cell constraints
  (default: `0.25 * max(pdims)`)
- `cc_rho`: Relative distance function for cell constraint spacing
  (default: `x -> ones(size(x, 1))`)
- `prot_layer`: Whether to add protection layers around cell constraints
  (default: false)
- `prot_d`: Protection layer distance functions (default: proportional to cell size)
- `face_constraints`: Vector of n×2 matrices defining fault/face constraint paths
  (default: empty)
- `fc_factor`: Relative grid size for face constraints (default: 1.0)
- `interpolate_fc`: Whether to interpolate face constraints (default: false)
- `circle_factor`: Circle size ratio for fault sites (default: 0.6, valid: (0.5, 1.0))
- `fc_rho`: Relative distance function for face constraint spacing
  (default: `x -> ones(size(x, 1))`)
- `fc_refinement`: Enable refinement around face constraints (default: false)
- `fc_eps`: Refinement transition distance around face constraints
  (default: `0.25 * max(pdims)`)
- `poly_bdr`: k×2 polygon boundary vertices (default: rectangular from `pdims`).
  Must have at least 3 vertices if specified.
- `suf_fc_cond`: Enforce the sufficient face constraint condition (default: true).
  When false, a less strict condition is used that removes reservoir sites
  closer than the constraint grid size.
- `max_iter`: Maximum DistMesh iterations (default: 500)
- `triangulation`: A [`TriangulationBackend`](@ref) instance for Delaunay
  triangulation (default: [`default_triangulation_backend()`](@ref))

# Returns
- `G::UnstructuredMesh`: Valid mesh with tagged cells and faces.
  `G.cells.tag` is `true` for cell constraint cells;
  `G.faces.tag` is `true` for face constraint faces.
- `pts`: Array of all Voronoi sites used.
- `F`: `SurfaceSitesResult` from face constraint processing.

# Example
```julia
fl = [Float64[0.2 0.2; 0.8 0.8]]
wl = [Float64[0.2 0.8; 0.8 0.2]]
G, pts, F = voronoi_mesh_2d(0.1, [1.0, 1.0];
    cell_constraints=wl, face_constraints=fl)
```

See also: [`composite_voronoi_mesh_2d`](@ref), [`clipped_voronoi_2d`](@ref)
"""
function voronoi_mesh_2d(
    res_grid_size::Float64,
    pdims::Vector{Float64};
    cell_constraints::Vector{Matrix{Float64}} = Matrix{Float64}[],
    cc_factor::Float64 = 1.0,
    interpolate_cc::Union{Bool,Vector{Bool}} = false,
    cc_refinement::Bool = false,
    cc_eps::Float64 = -1.0,
    cc_rho::Function = x -> ones(size(x, 1)),
    prot_layer::Bool = false,
    prot_d::Vector = Function[p -> ones(size(p, 1)) * res_grid_size / 10],
    face_constraints::Vector{Matrix{Float64}} = Matrix{Float64}[],
    fc_factor::Float64 = 1.0,
    interpolate_fc::Union{Bool,Vector{Bool}} = false,
    circle_factor::Float64 = 0.6,
    fc_rho::Function = x -> ones(size(x, 1)),
    fc_refinement::Bool = false,
    fc_eps::Float64 = -1.0,
    poly_bdr::Matrix{Float64} = zeros(0, 2),
    suf_fc_cond::Bool = true,
    max_iter::Int = 500,
    triangulation::TriangulationBackend = default_triangulation_backend()
)
    # ── Validate input ──
    @assert res_grid_size > 0
    @assert length(pdims) == 2
    @assert all(pdims .> 0)
    @assert 0.5 < circle_factor < 1.0

    # ── Compute grid sizes ──
    cc_grid_size = res_grid_size * cc_factor
    fc_grid_size = res_grid_size * fc_factor
    cc_rho_scaled = x -> cc_grid_size * cc_rho(x)
    @assert cc_grid_size > 0
    @assert fc_grid_size > 0

    # ── Handle polygon boundary ──
    if size(poly_bdr, 1) >= 3
        pdims = vec(maximum(poly_bdr, dims=1) .- minimum(poly_bdr, dims=1))
    elseif size(poly_bdr, 1) in (1, 2)
        error("Polygon must have at least 3 edges.")
    end

    # Default refinement transition distances
    if cc_eps < 0
        cc_eps = 0.25 * maximum(pdims)
    end
    if fc_eps < 0
        fc_eps = 0.25 * maximum(pdims)
    end

    # ── Handle scalar interpolation flags ──
    if isa(interpolate_cc, Bool)
        interpolate_cc = fill(interpolate_cc, max(length(cell_constraints), 1))
    end
    @assert length(interpolate_cc) >= length(cell_constraints)

    if isa(interpolate_fc, Bool)
        interpolate_fc = fill(interpolate_fc, max(length(face_constraints), 1))
    end
    @assert length(interpolate_fc) >= length(face_constraints)

    # Handle single protD
    if length(prot_d) == 1 && length(cell_constraints) > 1
        prot_d = fill(prot_d[1], length(cell_constraints))
    end
    if !isempty(cell_constraints)
        @assert length(prot_d) >= length(cell_constraints)
    end

    # ── Split constraints at intersections ──
    if !isempty(face_constraints) || !isempty(cell_constraints)
        fc_split, f_cut, fc_cut, fc_IC = split_at_intersections_2d(
            face_constraints, cell_constraints)
        cc_split, c_cut, cf_cut, cc_IC = split_at_intersections_2d(
            cell_constraints, face_constraints)

        # Map interpolation flags and protD
        if !isempty(fc_IC) && !isempty(face_constraints)
            interp_fl = [interpolate_fc[min(i, length(interpolate_fc))] for i in fc_IC]
        else
            interp_fl = Bool[]
            fc_split = face_constraints
            f_cut = zeros(Int, length(face_constraints))
            fc_cut = zeros(Int, length(face_constraints))
        end
        if !isempty(cc_IC) && !isempty(cell_constraints)
            interp_wp = [interpolate_cc[min(i, length(interpolate_cc))] for i in cc_IC]
            prot_d_mapped = [prot_d[min(i, length(prot_d))] for i in cc_IC]
        else
            interp_wp = Bool[]
            cc_split = cell_constraints
            c_cut = zeros(Int, length(cell_constraints))
            cf_cut = zeros(Int, length(cell_constraints))
            prot_d_mapped = prot_d
        end
    else
        fc_split = face_constraints
        cc_split = cell_constraints
        f_cut = zeros(Int, length(face_constraints))
        fc_cut = zeros(Int, length(face_constraints))
        c_cut = zeros(Int, length(cell_constraints))
        cf_cut = zeros(Int, length(cell_constraints))
        interp_fl = fill(false, length(face_constraints))
        interp_wp = fill(false, length(cell_constraints))
        prot_d_mapped = prot_d
    end

    # ── Handle point constraints (single-point cell constraints) ──
    for (idx, cc) in enumerate(cell_constraints)
        if size(cc, 1) == 1
            push!(cc_split, cc)
            push!(c_cut, 0)
            push!(cf_cut, 0)
            push!(interp_wp, false)
            if !isempty(prot_d)
                push!(prot_d_mapped, prot_d[min(idx, length(prot_d))])
            end
        end
    end

    # ── Calculate fault offset for se_ptn ──
    if cc_refinement
        f_len = cc_grid_size * fc_factor * 1.2
    else
        f_len = fc_grid_size
    end
    bisect_pnt = (f_len^2 - (circle_factor * f_len)^2 +
                  (circle_factor * f_len)^2) / (2 * f_len)
    fault_offset = sqrt(max(0, (circle_factor * f_len)^2 - bisect_pnt^2))

    # ── Create cell constraint sites ──
    well_pts = zeros(0, 2)
    w_gs = Float64[]
    prot_pts = zeros(0, 2)
    p_gs = Float64[]

    if !isempty(cc_split)
        se_ptn = zeros(length(cc_split), 2)
        for i in 1:length(cc_split)
            cf_i = cf_cut[min(i, length(cf_cut))]
            se_ptn[i, 1] = ((cf_i == 2 || cf_i == 3) ? 1.0 : 0.0) *
                           (1.0 + fault_offset / cc_grid_size)
            se_ptn[i, 2] = ((cf_i == 1 || cf_i == 3) ? 1.0 : 0.0) *
                           (1.0 + fault_offset / cc_grid_size)
        end

        well_pts, w_gs, prot_pts, p_gs = line_sites_2d(
            cc_split, cc_grid_size;
            se_ptn = se_ptn,
            cf_cut = cf_cut[1:min(length(cf_cut), length(cc_split))],
            c_cut = c_cut[1:min(length(c_cut), length(cc_split))],
            prot_layer = prot_layer,
            prot_d = prot_d_mapped[1:min(length(prot_d_mapped), length(cc_split))],
            cc_rho = cc_rho_scaled,
            interpolate_cc = interp_wp[1:min(length(interp_wp), length(cc_split))]
        )
    end

    # ── Create density functions ──
    if cc_refinement && size(well_pts, 1) > 0
        hresw = function (x)
            md = _min_pdist(x, well_pts)
            return min.(ones(size(x, 1)) ./ cc_factor, 1.2 .* exp.(md ./ cc_eps))
        end
        hfault = function (x)
            return cc_grid_size .* fc_factor .* hresw(x) .* fc_rho(x)
        end
    else
        hresw = x -> ones(size(x, 1)) ./ cc_factor
        hfault = x -> fc_grid_size .* fc_rho(x)
    end

    # ── Create surface (face constraint) sites ──
    F = if !isempty(fc_split)
        surface_sites_2d(
            fc_split, fc_grid_size;
            circle_factor = circle_factor,
            f_cut = f_cut[1:min(length(f_cut), length(fc_split))],
            fc_cut = fc_cut[1:min(length(fc_cut), length(fc_split))],
            interpolate_fc = interp_fl[1:min(length(interp_fl), length(fc_split))],
            dist_fun = hfault
        )
    else
        SurfaceSitesResult(
            zeros(0, 2), Float64[], Int[], [1],
            zeros(0, 2), Float64[], Int[], [1],
            Int[], [1], 0, zeros(0, 2)
        )
    end

    if fc_refinement && size(F.f_pts, 1) > 0
        hresf = function (x)
            md = _min_pdist(x, F.f_pts)
            return min.(ones(size(x, 1)) ./ fc_factor, 1.2 .* exp.(md ./ fc_eps))
        end
    else
        hresf = x -> ones(size(x, 1)) ./ cc_factor
    end

    # ── Set up domain ──
    if size(poly_bdr, 1) == 0
        x_max, y_max = pdims[1], pdims[2]
        rectangle = [0.0 0.0; x_max y_max]
        fd = p -> _drectangle(p, 0.0, x_max, 0.0, y_max)
        corners = [0.0 0.0; 0.0 y_max; x_max 0.0; x_max y_max]
        poly_bdr = [0.0 0.0; x_max 0.0; x_max y_max; 0.0 y_max]
    else
        rectangle = [minimum(poly_bdr, dims=1); maximum(poly_bdr, dims=1)]
        corners = copy(poly_bdr)
        fd = p -> _dpoly(p, poly_bdr)
    end

    # ── Remove tip sites outside domain ──
    if size(F.t_pts, 1) > 0
        inside = points_in_polygon(F.t_pts, poly_bdr)
        F.t_pts = F.t_pts[inside, :]
    end

    # ── Set up DistMesh density function ──
    if fc_refinement && cc_refinement
        ds = min(cc_grid_size, fc_grid_size)
        hres = x -> min.(hresf(x), hresw(x))
    elseif fc_refinement
        ds = fc_grid_size
        hres = hresf
    else
        ds = cc_grid_size
        hres = hresw
    end

    # ── Assemble fixed points and run DistMesh ──
    fixed_pts = vcat(F.f_pts, well_pts, prot_pts, F.t_pts, corners)

    pts_out, _ = distmesh_2d(fd, hres, ds, rectangle, fixed_pts;
        max_iter = max_iter, triangulation = triangulation)

    # ── Separate point types ──
    nf = size(F.f_pts, 1)
    nw = size(well_pts, 1)
    np = size(prot_pts, 1)
    nt = size(F.t_pts, 1)
    nc = size(corners, 1)
    n_fixed = nf + nw + np + nt + nc

    f_pts_out = pts_out[1:nf, :]
    w_pts_out = pts_out[nf+1:nf+nw, :]
    res_pts = pts_out[n_fixed+1:end, :]

    # ── Apply face constraint condition ──
    if suf_fc_cond && size(F.f_pts, 1) > 0
        res_pts, _ = surface_suf_cond_2d(res_pts, F)
    elseif size(F.f_pts, 1) > 0
        res_pts, _ = remove_conflict_points(res_pts, F.f_pts, F.f_Gs)
    end

    # ── Combine all points: fault, well, reservoir ──
    pts = vcat(f_pts_out, w_pts_out, res_pts)

    # ── Create the clipped Voronoi mesh ──
    G = clipped_voronoi_2d(pts, poly_bdr; triangulation = triangulation)

    # ── Tag face constraint faces ──
    G.faces.tag = fill(false, G.faces.num)
    if nf > 0 && length(F.f_c) >= nf
        for fi in 1:G.faces.num
            n1, n2 = G.faces.neighbors[fi, 1], G.faces.neighbors[fi, 2]
            if n1 == 0 || n2 == 0
                continue
            end
            # Both neighbors must be fault site cells
            if n1 <= nf && n2 <= nf
                c1 = F.f_c[n1]
                c2 = F.f_c[n2]
                if c1 == c2
                    G.faces.tag[fi] = true
                end
            end
        end
    end

    # ── Tag cell constraint cells ──
    G.cells.tag = fill(false, G.cells.num)
    if nw > 0
        well_start = nf + 1
        well_end = nf + nw
        for i in well_start:min(well_end, G.cells.num)
            G.cells.tag[i] = true
        end
    end

    return G, pts, F
end
