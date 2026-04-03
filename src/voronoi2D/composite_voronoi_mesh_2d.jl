"""
Construct a 2D composite Voronoi mesh with cell and face constraints.
Ported from MRST's compositePebiGrid2D.

This is the main entry point for creating conforming Voronoi meshes.
"""

"""
    composite_voronoi_mesh_2d(celldim, pdims; kwargs...)

Construct a 2D composite Voronoi mesh: a Cartesian background mesh refined
around face constraints (faults) and cell constraints (wells).

# Arguments
- `celldim`: [dx, dy] size of reservoir mesh cells
- `pdims`: [xmax, ymax] physical dimensions of the domain

# Keyword Arguments
- `cell_constraints`: Vector of n×2 matrices defining well/cell constraint paths
- `cc_factor`: Relative mesh size for cell constraints (default: 1.0)
- `interpolate_cc`: Whether to interpolate cell constraints (default: false)
- `mlqt_max_level`: Number of MLQT refinement levels (default: 0)
- `mlqt_level_steps`: Distance tolerance per level (default: -1)
- `cc_rho`: Relative distance function for cell constraint spacing
- `prot_layer`: Whether to add protection layers (default: false)
- `prot_d`: Protection layer distance functions
- `face_constraints`: Vector of n×2 matrices defining fault/face constraint paths
- `fc_factor`: Relative mesh size for face constraints (default: 1.0)
- `interpolate_fc`: Whether to interpolate face constraints (default: false)
- `circle_factor`: Circle size ratio for fault sites (default: 0.6, valid: (0.5, 1.0))
- `poly_bdr`: k×2 polygon boundary vertices (default: rectangular from pdims)
- `triangulation`: A [`TriangulationBackend`](@ref) instance for Delaunay triangulation
  (default: [`default_triangulation_backend()`](@ref))

# Returns
- `G::UnstructuredMesh`: Valid mesh with tagged cells and faces
- `pts`: All Voronoi sites used
- `F`: SurfaceSitesResult from face constraint processing

# Example
```julia
fl = [Float64[0.2 0.2; 0.8 0.8]]
wl = [Float64[0.2 0.8; 0.8 0.2]]
G, pts, F = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=wl, face_constraints=fl)
```
"""
function composite_voronoi_mesh_2d(
    celldim::Vector{Float64},
    pdims::Vector{Float64};
    cell_constraints::Vector{Matrix{Float64}} = Matrix{Float64}[],
    cc_factor::Float64 = 1.0,
    interpolate_cc::Union{Bool,Vector{Bool}} = false,
    mlqt_max_level::Int = 0,
    mlqt_level_steps::Float64 = -1.0,
    cc_rho::Function = x -> ones(size(x, 1)),
    prot_layer::Bool = false,
    prot_d::Vector = Function[p -> ones(size(p, 1)) * norm(celldim) / 10],
    face_constraints::Vector{Matrix{Float64}} = Matrix{Float64}[],
    fc_factor::Float64 = 1.0,
    interpolate_fc::Union{Bool,Vector{Bool}} = false,
    circle_factor::Float64 = 0.6,
    poly_bdr::Matrix{Float64} = zeros(0, 2),
    triangulation::TriangulationBackend = default_triangulation_backend()
)
    # Validate input
    @assert length(pdims) == 2
    @assert all(pdims .> 0)
    @assert all(celldim .> 0)
    @assert length(celldim) == 2
    @assert 0.5 < circle_factor < 1.0

    # Set mesh sizes
    cc_mesh_size = minimum(celldim) * cc_factor
    fc_mesh_size = minimum(celldim) * fc_factor

    # Handle scalar interpolate flags
    if isa(interpolate_cc, Bool)
        interpolate_cc = fill(interpolate_cc, max(length(cell_constraints), 1))
    end
    if isa(interpolate_fc, Bool)
        interpolate_fc = fill(interpolate_fc, max(length(face_constraints), 1))
    end

    # Handle single protD
    if length(prot_d) == 1 && length(cell_constraints) > 1
        prot_d = fill(prot_d[1], length(cell_constraints))
    end

    # Split face constraints and cell constraints at intersections
    if !isempty(face_constraints) || !isempty(cell_constraints)
        fc_split, f_cut, fc_cut, fc_IC = split_at_intersections_2d(face_constraints, cell_constraints)
        cc_split, c_cut, cf_cut, cc_IC = split_at_intersections_2d(cell_constraints, face_constraints)

        # Map interpolation flags
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

    # Handle point constraints
    for cc in cell_constraints
        if size(cc, 1) == 1
            push!(cc_split, cc)
            push!(c_cut, 0)
            push!(cf_cut, 0)
            push!(interp_wp, false)
            if !isempty(prot_d)
                push!(prot_d_mapped, prot_d[min(1, length(prot_d))])
            end
        end
    end

    # Calculate fault offset for se_ptn
    # bisect_pnt = (d² - R2² + R1²) / (2d) where d = R1 = R2 = fc_mesh_size terms
    # simplifies to fc_mesh_size / 2
    bisect_pnt = fc_mesh_size / 2
    fault_offset = sqrt(max(0, (circle_factor * fc_mesh_size)^2 - bisect_pnt^2))

    # Create cell constraint sites
    cc_pts = zeros(0, 2)
    c_gs = Float64[]
    prot_pts = zeros(0, 2)
    p_gs = Float64[]

    if !isempty(cc_split)
        se_ptn = zeros(length(cc_split), 2)
        for i in 1:length(cc_split)
            se_ptn[i, 1] = ((cf_cut[min(i, length(cf_cut))] == 2 ||
                             cf_cut[min(i, length(cf_cut))] == 3) ? 1.0 : 0.0) *
                           (1.0 + fault_offset / cc_mesh_size)
            se_ptn[i, 2] = ((cf_cut[min(i, length(cf_cut))] == 1 ||
                             cf_cut[min(i, length(cf_cut))] == 3) ? 1.0 : 0.0) *
                           (1.0 + fault_offset / cc_mesh_size)
        end

        cc_rho_scaled = x -> cc_mesh_size * cc_rho(x)

        cc_pts, c_gs, prot_pts, p_gs = line_sites_2d(
            cc_split, cc_mesh_size;
            se_ptn=se_ptn,
            cf_cut=cf_cut[1:min(length(cf_cut), length(cc_split))],
            c_cut=c_cut[1:min(length(c_cut), length(cc_split))],
            prot_layer=prot_layer,
            prot_d=prot_d_mapped[1:min(length(prot_d_mapped), length(cc_split))],
            cc_rho=cc_rho_scaled,
            interpolate_cc=interp_wp[1:min(length(interp_wp), length(cc_split))]
        )
    end

    # Create fault sites
    F = if !isempty(fc_split)
        surface_sites_2d(
            fc_split, fc_mesh_size;
            circle_factor=circle_factor,
            f_cut=f_cut[1:min(length(f_cut), length(fc_split))],
            fc_cut=fc_cut[1:min(length(fc_cut), length(fc_split))],
            interpolate_fc=interp_fl[1:min(length(interp_fl), length(fc_split))],
            dist_fun=x -> fc_mesh_size * ones(size(x, 1))
        )
    else
        SurfaceSitesResult(
            zeros(0, 2), Float64[], Int[], [1],
            zeros(0, 2), Float64[], Int[], [1],
            Int[], [1], 0, zeros(0, 2)
        )
    end

    # Create reservoir background mesh
    if size(poly_bdr, 1) == 0
        dx = pdims[1] / ceil(pdims[1] / celldim[1])
        dy = pdims[2] / ceil(pdims[2] / celldim[2])
        vx = collect(0:dx:pdims[1])
        vy = collect(0:dy:pdims[2])
        poly_bdr = [0.0 0.0; pdims[1] 0.0; pdims[1] pdims[2]; 0.0 pdims[2]]
    elseif size(poly_bdr, 1) < 3
        error("Polygon must have at least 3 edges")
    else
        l_dim_min = minimum(poly_bdr, dims=1)
        l_dim_max = maximum(poly_bdr, dims=1)
        dx = (l_dim_max[1] - l_dim_min[1]) / ceil((l_dim_max[1] - l_dim_min[1]) / celldim[1])
        dy = (l_dim_max[2] - l_dim_min[2]) / ceil((l_dim_max[2] - l_dim_min[2]) / celldim[2])
        vx = collect(l_dim_min[1]:dx:l_dim_max[1])
        vy = collect(l_dim_min[2]:dy:l_dim_max[2])
    end

    # Create meshgrid
    res_pts_init = zeros(length(vx) * length(vy), 2)
    idx = 1
    for y in vy, x in vx
        res_pts_init[idx, :] = [x, y]
        idx += 1
    end

    # Remove points outside polygon (if not rectangular)
    if size(poly_bdr, 1) >= 3 && !(size(poly_bdr, 1) == 4 &&
                                     poly_bdr[1, :] ≈ [0.0, 0.0] &&
                                     poly_bdr[3, :] ≈ pdims)
        inside = points_in_polygon(res_pts_init, poly_bdr)
        res_pts_init = res_pts_init[inside, :]
    end

    # Remove tip sites outside domain
    if size(F.t_pts, 1) > 0
        inside = points_in_polygon(F.t_pts, poly_bdr)
        F.t_pts = F.t_pts[inside, :]
    end

    # MLQT refinement
    res_pts = if !isempty(cc_pts) && mlqt_max_level > 0
        refined = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
        for i in 1:size(res_pts_init, 1)
            dt = mlqt_level_steps > 0 ? mlqt_level_steps : -1.0
            sub = mlqt(
                res_pts_init[i, :], cc_pts, celldim;
                level=1, max_level=mlqt_max_level, dist_tol=dt
            )
            append!(refined, sub)
        end
        reduce(vcat, [r[1]' for r in refined])
    else
        res_pts_init
    end

    # Remove conflict points
    if size(cc_pts, 1) > 0
        res_pts, _ = remove_conflict_points(res_pts, cc_pts, c_gs)
    end
    if size(prot_pts, 1) > 0
        res_pts, _ = remove_conflict_points(res_pts, prot_pts, p_gs)
    end
    if size(F.f_pts, 1) > 0
        res_pts, _ = remove_conflict_points(res_pts, F.f_pts, F.f_Gs)
    end
    if size(F.c_CC, 1) > 0
        res_pts, _ = remove_conflict_points(res_pts, F.c_CC, F.c_R)
    end

    # Combine all points
    pts = vcat(F.f_pts, cc_pts, prot_pts, F.t_pts, res_pts)

    # Create mesh
    G = clipped_voronoi_2d(pts, poly_bdr; triangulation=triangulation)

    # Tag fault faces
    G.faces.tag = fill(false, G.faces.num)

    # Tag constrained cells
    G.cells.tag = fill(false, G.cells.num)
    if size(cc_pts, 1) > 0
        n_fault_pts = size(F.f_pts, 1)
        well_start = n_fault_pts + 1
        well_end = n_fault_pts + size(cc_pts, 1)
        for i in well_start:min(well_end, G.cells.num)
            G.cells.tag[i] = true
        end
    end

    return G, pts, F
end
