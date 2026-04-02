"""
Multi-level quadtree refinement.
Ported from MRST's mlqt function.
"""

"""
    mlqt(cell_center, bndr, cell_size; level=1, max_level=0, dist_tol=-1.0)

Perform multi-level mesh refinement around boundary points.

# Arguments
- `cell_center`: Center coordinate of cell (1×2 or 1×3)
- `bndr`: n×d array of coordinates where cell should be refined
- `cell_size`: Size of the cell [dx, dy] or [dx, dy, dz]

# Keyword Arguments
- `level`: Current refinement level (default: 1)
- `max_level`: Maximum refinement level (default: 0)
- `dist_tol`: Distance tolerance for refinement; -1 uses default (default: -1)

# Returns
- Vector of tuples `(center, size)` for refined cells
"""
function mlqt(
    cell_center::Vector{Float64},
    bndr::Matrix{Float64},
    cell_size::Vector{Float64};
    level::Int = 1,
    max_level::Int = 0,
    dist_tol::Union{Float64,Vector{Float64}} = -1.0
)
    # Is recursion finished?
    if level > max_level
        return [(copy(cell_center), copy(cell_size))]
    end

    # Set distance tolerance
    if isa(dist_tol, Float64) || length(dist_tol) == 1
        dt = isa(dist_tol, Float64) ? dist_tol : dist_tol[1]
        if dt <= 0
            dist_tol_val = 2.5 * cell_size / 2
        else
            dist_tol_val = fill(dt, length(cell_size))
        end
        dist_next = dist_tol_val / 2
    else
        dist_tol_val = fill(dist_tol[level], length(cell_size))
        dist_next = dist_tol
    end

    # Should cell be refined?
    should_refine = false
    for i in 1:size(bndr, 1)
        if all(abs.(cell_center .- bndr[i, :]) .<= dist_tol_val)
            should_refine = true
            break
        end
    end

    if should_refine
        dx = cell_size / 4
        new_size = cell_size / 2

        # 2D quadtree shifts
        shifts = [
            dx[1] dx[2];
            dx[1] -dx[2];
            -dx[1] -dx[2];
            -dx[1] dx[2]
        ]

        if length(cell_center) == 3
            shifts_2d = copy(shifts)
            shifts = vcat(
                hcat(shifts_2d, fill(dx[3], 4)),
                hcat(shifts_2d, fill(-dx[3], 4))
            )
        end

        results = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
        for i in 1:size(shifts, 1)
            sub_results = mlqt(
                cell_center .+ shifts[i, :],
                bndr, new_size;
                level=level + 1,
                max_level=max_level,
                dist_tol=isa(dist_next, Vector) ? dist_next : [dist_next[1]]
            )
            append!(results, sub_results)
        end
        return results
    else
        return [(copy(cell_center), copy(cell_size))]
    end
end
