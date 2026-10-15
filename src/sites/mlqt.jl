# Multilevel quadtree/octree refinement of Cartesian background sites around
# cell constraints (UPR: mlqt).

"""
    mlqt(cell_center, constraint_sites, cell_size; max_level = 0, dist_tol = -1)

Refine the Cartesian cell with center `cell_center` and size `cell_size` (a
vector, one entry per dimension) into 4 (2D) or 8 (3D) children if any of the
`constraint_sites` is within `dist_tol` (per dimension) of the center, and
recurse up to `max_level` times.

`dist_tol <= 0` gives the default `1.25*cell_size`, halved for each level. A
vector of length `max_level` gives the tolerance per level.

Returns the vector of `(center, size)` of the leaf cells.
"""
function mlqt(cell_center::SVector{D, Float64}, constraint_sites::AbstractVector{SVector{D, Float64}}, cell_size;
        max_level = 0, dist_tol = -1, level = 1, tree = nothing) where {D}
    cell_size = SVector{D, Float64}(cell_size)
    all(cell_size .> 0) || throw(ArgumentError("cell_size must be positive"))
    if level > max_level || isempty(constraint_sites)
        return [(cell_center, cell_size)]
    end
    if isnothing(tree)
        tree = KDTree(collect(constraint_sites))
    end
    if dist_tol isa AbstractVector && length(dist_tol) > 1
        tol = fill(Float64(dist_tol[level]), D)
        next_tol = dist_tol
    else
        t = dist_tol isa AbstractVector ? only(dist_tol) : dist_tol
        # The default tolerance follows the cell size at each level
        tol = t <= 0 ? 2.5 .* cell_size ./ 2 : fill(Float64(t), D)
        next_tol = t <= 0 ? -1 : t / 2
    end
    tol = SVector{D, Float64}(tol)
    refine = false
    for i in inrange(tree, cell_center, norm(tol))
        if all(abs.(constraint_sites[i] - cell_center) .<= tol)
            refine = true
            break
        end
    end
    refine || return [(cell_center, cell_size)]
    dx = cell_size / 4
    out = Tuple{SVector{D, Float64}, SVector{D, Float64}}[]
    for signs in Iterators.product(ntuple(_ -> (1, -1), D)...)
        shift = SVector{D, Float64}(signs) .* dx
        append!(out, mlqt(cell_center + shift, constraint_sites, cell_size / 2;
            max_level = max_level, dist_tol = next_tol, level = level + 1, tree = tree))
    end
    return out
end
