"""
Place Voronoi sites along cell constraints (wells).
Ported from MRST's lineSites2D.
"""

"""
    line_sites_2d(cell_constraints, cc_mesh_size; kwargs...)

Place sites along given lines (cell constraints / wells).

# Arguments
- `cell_constraints`: Vector of n×2 matrices defining constraint paths
- `cc_mesh_size`: Desired distance between sites along constraints

# Keyword Arguments
- `se_ptn`: Start/end position offsets per constraint (n×2 matrix)
- `cf_cut`: Cut indicators for constraint-fault intersections
- `c_cut`: Cut indicators for constraint-constraint intersections
- `prot_layer`: Whether to add protection layers (default: false)
- `prot_d`: Protection layer distance functions
- `cc_rho`: Distance function for constraint spacing
- `interpolate_cc`: Whether to interpolate constraints (default: false)

# Returns
- `cc_pts`: All generated constraint sites (n×2)
- `c_gs`: Mesh spacing for each site
- `prot_pts`: Protection layer sites (n×2)
- `p_gs`: Protection layer spacings
"""
function line_sites_2d(
    cell_constraints::Vector{Matrix{Float64}},
    cc_mesh_size::Float64;
    se_ptn::Matrix{Float64} = zeros(length(cell_constraints), 2),
    cf_cut::Vector{Int} = zeros(Int, length(cell_constraints)),
    c_cut::Vector{Int} = zeros(Int, length(cell_constraints)),
    prot_layer::Bool = false,
    prot_d::Vector = Function[p -> ones(size(p, 1)) * cc_mesh_size / 10],
    cc_rho::Function = x -> cc_mesh_size * ones(size(x, 1)),
    interpolate_cc::Union{Bool,Vector{Bool}} = false
)
    nc = length(cell_constraints)

    # Handle scalar interpolate_cc
    if isa(interpolate_cc, Bool)
        interpolate_cc = fill(interpolate_cc, nc)
    end

    # Handle single protD
    if length(prot_d) == 1 && nc > 1
        prot_d = fill(prot_d[1], nc)
    end

    cc_pts = zeros(0, 2)
    c_gs = Float64[]
    prot_pts = zeros(0, 2)
    p_gs = Float64[]

    for i in 1:nc
        constraint = cell_constraints[i]

        if size(constraint, 1) == 1
            # Point constraint
            p = copy(constraint)
            well_space = [cc_mesh_size]
        else
            p = inter_line_path(constraint, cc_rho, cc_mesh_size,
                se_ptn[i, :], interpolate_cc[i])
            if isempty(p) || size(p, 1) == 0
                continue
            end
            ds = vec(sqrt.(sum(diff(p, dims=1) .^ 2, dims=2)))
            well_space = min.([ds[1]; ds], [ds; ds[end]])
        end

        # Apply cfCut
        keep = collect(1:size(p, 1))
        if cf_cut[i] == 1
            keep = keep[1:end-1]
        elseif cf_cut[i] == 2
            keep = keep[2:end]
        elseif cf_cut[i] == 3
            keep = keep[2:end-1]
        end

        keep_prot = copy(keep)
        if c_cut[i] == 1
            keep_prot = keep_prot[1:end-1]
        elseif c_cut[i] == 2
            keep_prot = keep_prot[2:end]
        elseif c_cut[i] == 3
            keep_prot = keep_prot[2:end-1]
        end

        if !isempty(keep)
            append!(c_gs, well_space[keep])
            cc_pts = vcat(cc_pts, p[keep, :])
        end

        # Add protection layer
        if size(constraint, 1) > 1 && prot_layer && !isempty(keep_prot)
            if length(keep_prot) == 1
                pK = vcat(p[keep_prot, :], constraint[end:end, :])
            else
                pK = p[keep_prot, :]
            end

            # Calculate normals
            new_n = diff(pK, dims=1)
            new_n = hcat(new_n[:, 2], -new_n[:, 1])
            norms = sqrt.(sum(new_n .^ 2, dims=2))
            new_n = new_n ./ norms

            if length(keep_prot) > 1
                new_n = vcat(new_n, new_n[end:end, :])
            else
                pK = pK[1:1, :]
            end

            d_vals = prot_d[i](pK)
            d_mat = hcat(d_vals, d_vals)

            prot_pts = vcat(prot_pts, pK .+ new_n .* d_mat, pK .- new_n .* d_mat)
            append!(p_gs, vec(d_mat))
        end
    end

    # Remove duplicates
    if size(cc_pts, 1) > 0
        cc_pts_rounded = round.(cc_pts .* 1e13) ./ 1e13
        seen = Set{Tuple{Float64,Float64}}()
        keep_idx = Int[]
        for i in 1:size(cc_pts_rounded, 1)
            key = (cc_pts_rounded[i, 1], cc_pts_rounded[i, 2])
            if !(key in seen)
                push!(seen, key)
                push!(keep_idx, i)
            end
        end
        cc_pts = cc_pts[keep_idx, :]
        c_gs = c_gs[keep_idx]
    end

    return cc_pts, c_gs, prot_pts, p_gs
end
