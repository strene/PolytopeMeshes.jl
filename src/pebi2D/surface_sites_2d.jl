"""
Place Voronoi sites on both sides of face constraints (faults/fractures).
Ported from MRST's surfaceSites2D.
"""

"""
    SurfaceSitesResult

Result struct from `surface_sites_2d`.

# Fields
- `f_pts`: Face constraint site coordinates
- `f_Gs`: Grid spacing for each site
- `f_c`: Map from sites to circles
- `f_cPos`: Position pointers for site-to-circle map
- `c_CC`: Circle center coordinates
- `c_R`: Circle radii
- `c_f`: Map from circles to sites
- `c_fPos`: Position pointers for circle-to-site map
- `c_l`: Map from circles to face constraints
- `l_fPos`: Map from face constraints to sites (position pointers)
- `l_nFault`: Number of face constraints
- `t_pts`: Tip site coordinates
"""
mutable struct SurfaceSitesResult
    f_pts::Matrix{Float64}
    f_Gs::Vector{Float64}
    f_c::Vector{Int}
    f_cPos::Vector{Int}
    c_CC::Matrix{Float64}
    c_R::Vector{Float64}
    c_f::Vector{Int}
    c_fPos::Vector{Int}
    c_l::Vector{Int}
    l_fPos::Vector{Int}
    l_nFault::Int
    t_pts::Matrix{Float64}
end

"""
    surface_sites_2d(face_constraints, fc_grid_size; kwargs...)

Place surface sites on both sides of face constraints (faults).

# Arguments
- `face_constraints`: Vector of n×2 matrices defining fault paths
- `fc_grid_size`: Desired distance between fault sites

# Keyword Arguments
- `circle_factor`: Ratio of circle radius to grid size (default: 0.6, valid: (0.5, 1.0))
- `f_cut`: L1-L1 intersection cut indicators
- `fc_cut`: L1-L2 intersection cut indicators
- `interpolate_fc`: Whether to interpolate constraints
- `dist_fun`: Distance function for adaptive spacing

# Returns
- `F::SurfaceSitesResult`: Struct containing all site information
"""
function surface_sites_2d(
    face_constraints::Vector{Matrix{Float64}},
    fc_grid_size::Float64;
    circle_factor::Float64 = 0.6,
    f_cut::Vector{Int} = zeros(Int, length(face_constraints)),
    fc_cut::Vector{Int} = zeros(Int, length(face_constraints)),
    interpolate_fc::Union{Bool,Vector{Bool}} = false,
    dist_fun::Function = x -> fc_grid_size * ones(size(x, 1))
)
    n_faults = length(face_constraints)

    if isa(interpolate_fc, Bool)
        interpolate_fc = fill(interpolate_fc, n_faults)
    end

    # Initialize result
    F = SurfaceSitesResult(
        zeros(0, 2),  # f_pts
        Float64[],    # f_Gs
        Int[],        # f_c
        [1],          # f_cPos
        zeros(0, 2),  # c_CC
        Float64[],    # c_R
        Int[],        # c_f
        [1],          # c_fPos
        Int[],        # c_l
        [1],          # l_fPos
        n_faults,
        zeros(0, 2)   # t_pts
    )

    for i in 1:n_faults
        fc = face_constraints[i]
        se_ptn = [
            0.5 * ((fc_cut[i] == 2 || fc_cut[i] == 3) ? 1.0 : 0.0),
            0.5 * ((fc_cut[i] == 1 || fc_cut[i] == 3) ? 1.0 : 0.0)
        ]

        pts, grid_spacing, circ_center, circ_radius =
            _face_constraint_pts(fc, fc_grid_size, circle_factor,
                f_cut[i], se_ptn, dist_fun, interpolate_fc[i])

        nl = div(size(pts, 1), 2)
        if nl == 0
            push!(F.l_fPos, F.l_fPos[end])
            continue
        end

        append!(F.f_Gs, grid_spacing)
        push!(F.l_fPos, size(F.f_pts, 1) + 1 + size(pts, 1))
        F.f_pts = vcat(F.f_pts, pts)
        F.c_CC = vcat(F.c_CC, circ_center)
        append!(F.c_R, circ_radius)

        # Build simple circle-site mappings
        nc = size(circ_center, 1)
        n_existing = length(F.c_fPos) - 1

        for j in 1:nc
            # Each circle creates 2 sites (left and right)
            push!(F.c_f, size(F.f_pts, 1) - 2 * nl + j)
            push!(F.c_f, size(F.f_pts, 1) - nl + j)
            push!(F.c_fPos, F.c_fPos[end] + 2)
            push!(F.c_l, i)
        end

        # Site to circle mapping
        for j in 1:(2 * nl)
            circ_idx = min(div(j - 1, 1) + 1, nc)
            if j <= nl
                circ_idx = j
            else
                circ_idx = j - nl
            end
            push!(F.f_c, n_existing + circ_idx)
            push!(F.f_cPos, F.f_cPos[end] + 1)
        end
    end

    # Add tip sites
    if size(F.c_CC, 1) > 0
        _add_tip_sites!(F, face_constraints, fc_cut, f_cut)
    end

    return F
end

"""
Generate face constraint points for a single fault.
Returns: pts, grid_spacing, circ_center, circ_radius
"""
function _face_constraint_pts(
    fc::Matrix{Float64},
    frac_ds::Float64,
    circle_factor::Float64,
    is_cut::Int,
    se_ptn::Vector{Float64},
    fh::Function,
    interp_fl::Bool
)
    @assert 0.5 < circle_factor < 1.0
    @assert size(fc, 1) > 1 && size(fc, 2) == 2

    # Interpolate FC line to get circle centers
    circ_center = inter_line_path(fc, fh, frac_ds, se_ptn, interp_fl)
    num_frac_pts = size(circ_center, 1) - 1

    # Test if face constraint is too short
    if num_frac_pts == 1
        d = norm(circ_center[2, :] - circ_center[1, :])
        mid = (circ_center[2, :] + circ_center[1, :]) / 2
        if d < 0.8 * fh(mid')[1]
            return zeros(0, 2), Float64[], zeros(0, 2), Float64[]
        end
    end

    if size(circ_center, 1) < 2
        return zeros(0, 2), Float64[], zeros(0, 2), Float64[]
    end

    # Calculate line lengths and circle radii
    line_length = vec(sqrt.(sum((circ_center[2:end, :] .- circ_center[1:end-1, :]) .^ 2, dims=2)))

    circ_radius = circle_factor .* vcat(
        [line_length[1]],
        (line_length[1:end-1] .+ line_length[2:end]) ./ 2,
        [line_length[end]]
    )

    # Adjust radii at cut ends
    if is_cut == 1
        circ_radius[end] = fh(circ_center[end:end, :])[1] * circle_factor
    elseif is_cut == 2
        circ_radius[1] = fh(circ_center[1:1, :])[1] * circle_factor
    elseif is_cut == 3
        circ_radius[1] = fh(circ_center[1:1, :])[1] * circle_factor
        circ_radius[end] = fh(circ_center[end:end, :])[1] * circle_factor
    end

    # Calculate circle intersections (bisection points)
    bisect_pnt = (line_length .^ 2 .- circ_radius[2:end] .^ 2 .+ circ_radius[1:end-1] .^ 2) ./ (2 .* line_length)

    # Prevent NaN from negative values under sqrt
    inner = circ_radius[1:end-1] .^ 2 .- bisect_pnt .^ 2
    inner = max.(inner, 0.0)
    fault_offset = sqrt.(inner)

    # Unit vectors along and normal to fault
    n1 = (circ_center[2:end, :] .- circ_center[1:end-1, :]) ./ line_length
    n2 = hcat(-n1[:, 2], n1[:, 1])

    # Set FC sites on both sides
    left = circ_center[1:end-1, :] .+ bisect_pnt .* n1 .+ fault_offset .* n2
    right = circ_center[1:end-1, :] .+ bisect_pnt .* n1 .- fault_offset .* n2

    pts = vcat(right, left)
    grid_spacing = vcat(2 .* fault_offset, 2 .* fault_offset)

    return pts, grid_spacing, circ_center, circ_radius
end

"""
Add tip sites at the ends of face constraints.
"""
function _add_tip_sites!(F::SurfaceSitesResult,
    face_constraints::Vector{Matrix{Float64}},
    fc_cut::Vector{Int},
    f_cut::Vector{Int})

    tip_pts = zeros(0, 2)

    # Track cumulative circle index per fault
    circ_offset = 0

    for i in 1:F.l_nFault
        fc = face_constraints[i]
        nl = div(F.l_fPos[i + 1] - F.l_fPos[i], 2)
        if nl == 0
            continue
        end
        nc = nl + 1  # number of circles for this fault

        # Add tip at start (if no cut)
        tip_at_start = !((fc_cut[i] == 2 || fc_cut[i] == 3) ||
                         (f_cut[i] == 2 || f_cut[i] == 3))
        tip_at_end = !((fc_cut[i] == 1 || fc_cut[i] == 3) ||
                       (f_cut[i] == 1 || f_cut[i] == 3))

        if tip_at_start && nc >= 2
            c1 = circ_offset + 1
            c2 = circ_offset + 2
            t_vec = F.c_CC[c1, :] .- F.c_CC[c2, :]
            t_len = norm(t_vec)
            if t_len > 0
                t_vec ./= t_len
                tip_pts = vcat(tip_pts, (F.c_CC[c1, :] .+ t_vec .* F.c_R[c1])')
            end
        end

        if tip_at_end && nc >= 2
            c1 = circ_offset + nc - 1
            c2 = circ_offset + nc
            t_vec = F.c_CC[c2, :] .- F.c_CC[c1, :]
            t_len = norm(t_vec)
            if t_len > 0
                t_vec ./= t_len
                tip_pts = vcat(tip_pts, (F.c_CC[c2, :] .+ t_vec .* F.c_R[c2])')
            end
        end

        circ_offset += nc
    end

    F.t_pts = vcat(F.t_pts, tip_pts)
end
