"""
2D mesh generation using the DistMesh algorithm.

Ported from the DistMesh algorithm by Per-Olof Persson and Gilbert Strang,
as adapted in MRST's UPR module.

Reference:
  Persson, P.-O. and Strang, G. (2004).
  "A simple mesh generator in MATLAB."
  SIAM Review, 46(2), 329–345.
"""

"""
    _drectangle(p, x1, x2, y1, y2)

Signed distance function for a rectangle [x1,x2] × [y1,y2].
Returns negative values inside, positive outside.

# Arguments
- `p`: n×2 matrix of point coordinates
- `x1, x2`: x-range of rectangle
- `y1, y2`: y-range of rectangle
"""
function _drectangle(p::Matrix{Float64}, x1::Float64, x2::Float64,
    y1::Float64, y2::Float64)
    return -min.(
        min.(p[:, 2] .- y1, y2 .- p[:, 2]),
        min.(p[:, 1] .- x1, x2 .- p[:, 1])
    )
end

"""
    _dist_to_segment(px, py, ax, ay, bx, by)

Compute the distance from point (px,py) to the line segment from (ax,ay) to (bx,by).
"""
function _dist_to_segment(px::Float64, py::Float64,
    ax::Float64, ay::Float64, bx::Float64, by::Float64)
    dx = bx - ax
    dy = by - ay
    len_sq = dx * dx + dy * dy
    if len_sq < eps()
        return sqrt((px - ax)^2 + (py - ay)^2)
    end
    t = clamp(((px - ax) * dx + (py - ay) * dy) / len_sq, 0.0, 1.0)
    proj_x = ax + t * dx
    proj_y = ay + t * dy
    return sqrt((px - proj_x)^2 + (py - proj_y)^2)
end

"""
    _dpoly(p, pv)

Signed distance function for a polygon.
Returns negative values inside, positive outside.

# Arguments
- `p`: n×2 matrix of point coordinates
- `pv`: m×2 matrix of polygon vertices (not necessarily closed; closing is handled internally)
"""
function _dpoly(p::Matrix{Float64}, pv::Matrix{Float64})
    np = size(p, 1)
    nv = size(pv, 1)
    d = zeros(np)

    for i in 1:np
        px, py = p[i, 1], p[i, 2]
        min_dist = Inf

        for j in 1:nv
            j_next = mod(j, nv) + 1
            dist = _dist_to_segment(px, py,
                pv[j, 1], pv[j, 2], pv[j_next, 1], pv[j_next, 2])
            if dist < min_dist
                min_dist = dist
            end
        end

        # Sign: negative inside, positive outside
        inside = point_in_polygon(px, py, pv)
        d[i] = inside ? -min_dist : min_dist
    end

    return d
end

"""
    _min_pdist(x, y)

Compute the minimum Euclidean distance from each row of `x` to the nearest row of `y`.

# Arguments
- `x`: n×2 matrix of query points
- `y`: m×2 matrix of reference points

# Returns
- Vector of length n with minimum distances
"""
function _min_pdist(x::Matrix{Float64}, y::Matrix{Float64})
    nx = size(x, 1)
    ny = size(y, 1)
    d = zeros(nx)

    for i in 1:nx
        min_d_sq = Inf
        for j in 1:ny
            d_sq = (x[i, 1] - y[j, 1])^2 + (x[i, 2] - y[j, 2])^2
            if d_sq < min_d_sq
                min_d_sq = d_sq
            end
        end
        d[i] = sqrt(min_d_sq)
    end

    return d
end

"""
    distmesh_2d(fd, fh, h0, bbox, p_fix; kwargs...)

Generate a 2D point distribution using the DistMesh algorithm (force-based
equilibration with Delaunay retriangulation).

# Arguments
- `fd`: Signed distance function `fd(p::Matrix{Float64}) -> Vector{Float64}`.
  Returns negative values inside the domain, positive outside.
- `fh`: Element size function `fh(p::Matrix{Float64}) -> Vector{Float64}`.
  Returns the desired relative edge length at each point.
- `h0`: Initial (base) edge length for the background grid.
- `bbox`: Bounding box `[x_min y_min; x_max y_max]` (2×2 matrix).
- `p_fix`: Fixed points that must appear in the output (n×2 matrix).
  These points are never moved or removed. They appear first in the output.

# Keyword Arguments
- `max_iter`: Maximum number of iterations (default: 500)
- `triangulation`: A [`TriangulationBackend`](@ref) instance for Delaunay
  triangulation (default: [`default_triangulation_backend()`](@ref))

# Returns
- `p`: Final point coordinates (n×2 matrix). Fixed points are the first `nfix` rows.
- `t`: Final Delaunay triangulation (m×3 matrix of vertex indices).
"""
function distmesh_2d(fd::Function, fh::Function, h0::Float64,
    bbox::Matrix{Float64}, p_fix::Matrix{Float64};
    max_iter::Int = 500,
    triangulation::TriangulationBackend = default_triangulation_backend())

    # Standard DistMesh algorithm parameters (Persson & Strang, 2004).
    # These are well-established defaults and are not exposed as keyword arguments
    # to keep the API simple.
    dptol = 0.001       # convergence tolerance: max relative displacement
    ttol = 0.1          # retriangulation threshold: max relative displacement
    Fscale = 1.2        # force scaling: ratio of desired to actual edge length
    deltat = 0.2        # time step for point movement (damping factor)
    geps = 0.001 * h0   # tolerance for inside/outside tests
    deps = sqrt(eps()) * h0  # step size for numerical gradient

    nfix = size(p_fix, 1)

    # 1. Create initial hexagonal point distribution in bounding box
    x_vals = collect(bbox[1, 1]:h0:bbox[2, 1])
    y_step = h0 * sqrt(3) / 2
    y_vals = collect(bbox[1, 2]:y_step:bbox[2, 2])

    px = Float64[]
    py = Float64[]
    for (iy, y) in enumerate(y_vals)
        offset = iseven(iy) ? h0 / 2 : 0.0
        for x in x_vals
            push!(px, x + offset)
            push!(py, y)
        end
    end

    if isempty(px)
        return copy(p_fix), zeros(Int, 0, 3)
    end

    p = hcat(px, py)

    # 2. Remove points outside domain
    d = fd(p)
    p = p[vec(d) .< geps, :]

    if size(p, 1) == 0
        return copy(p_fix), zeros(Int, 0, 3)
    end

    # Apply probability rejection based on element size function
    r0 = 1.0 ./ vec(fh(p)) .^ 2
    max_r0 = maximum(r0)
    if max_r0 > 0
        p = p[rand(size(p, 1)) .< r0 ./ max_r0, :]
    end

    if size(p, 1) == 0
        return copy(p_fix), zeros(Int, 0, 3)
    end

    # Remove generated points that are too close to fixed points
    if nfix > 0
        keep = trues(size(p, 1))
        min_dist_sq = (h0 / 2)^2
        for i in 1:size(p, 1)
            for j in 1:nfix
                if (p[i, 1] - p_fix[j, 1])^2 + (p[i, 2] - p_fix[j, 2])^2 < min_dist_sq
                    keep[i] = false
                    break
                end
            end
        end
        p = p[keep, :]
    end

    # Prepend fixed points
    p = vcat(p_fix, p)
    N = size(p, 1)

    p_old = fill(Inf, N, 2)
    bars = zeros(Int, 0, 2)

    for iter in 1:max_iter
        N = size(p, 1)

        # 3. Retriangulation if any point moved significantly
        max_move = sqrt(maximum(sum((p .- p_old[1:min(N, size(p_old, 1)), :]) .^ 2, dims=2)))
        if max_move / h0 > ttol
            p_old = copy(p)

            t = triangulate_2d(triangulation, p)

            if size(t, 1) == 0
                break
            end

            # Remove exterior triangles
            pmid = (p[t[:, 1], :] .+ p[t[:, 2], :] .+ p[t[:, 3], :]) ./ 3
            d_mid = fd(pmid)
            t = t[vec(d_mid) .< -geps, :]

            if size(t, 1) == 0
                break
            end

            # Extract unique edges (bars)
            bars_all = vcat(t[:, [1, 2]], t[:, [1, 3]], t[:, [2, 3]])
            bars_sorted = sort(bars_all, dims=2)
            bars = unique(bars_sorted, dims=1)
        end

        if size(bars, 1) == 0
            continue
        end

        # 4. Compute bar forces
        barvec = p[bars[:, 1], :] .- p[bars[:, 2], :]
        L = vec(sqrt.(sum(barvec .^ 2, dims=2)))

        hbars = vec(fh((p[bars[:, 1], :] .+ p[bars[:, 2], :]) ./ 2))
        sum_L2 = sum(L .^ 2)
        sum_h2 = sum(hbars .^ 2)
        if sum_h2 < eps()
            break
        end
        L0 = hbars .* Fscale .* sqrt(sum_L2 / sum_h2)

        # Spring force: repulsive only (zero if bar is shorter than desired)
        F = max.(L0 .- L, 0.0)

        # Avoid division by zero for zero-length bars
        L_safe = max.(L, eps())
        Fvec = (F ./ L_safe) .* barvec

        # 5. Sum forces at each node
        Ftot = zeros(N, 2)
        for k in 1:size(bars, 1)
            i, j = bars[k, 1], bars[k, 2]
            Ftot[i, 1] += Fvec[k, 1]
            Ftot[i, 2] += Fvec[k, 2]
            Ftot[j, 1] -= Fvec[k, 1]
            Ftot[j, 2] -= Fvec[k, 2]
        end

        # Zero force at fixed points
        Ftot[1:nfix, :] .= 0.0

        # Update positions
        p .+= deltat .* Ftot

        # 6. Project points outside domain back to boundary
        d = vec(fd(p))
        ix = findall(i -> i > nfix && d[i] > 0, 1:N)

        if !isempty(ix)
            p_ix = p[ix, :]
            d_ix = d[ix]

            # Numerical gradient of distance function
            px_plus = copy(p_ix)
            px_plus[:, 1] .+= deps
            py_plus = copy(p_ix)
            py_plus[:, 2] .+= deps

            dg_x = (vec(fd(px_plus)) .- d_ix) ./ deps
            dg_y = (vec(fd(py_plus)) .- d_ix) ./ deps

            # Project back: p -= d * grad(fd)
            p[ix, 1] .-= d_ix .* dg_x
            p[ix, 2] .-= d_ix .* dg_y
        end

        # 7. Check convergence
        d_check = vec(fd(p))
        inner = d_check .< -geps
        if !any(inner)
            break
        end
        inner_disp = deltat .* Ftot[inner, :]
        if size(inner_disp, 1) > 0
            max_disp = maximum(sqrt.(sum(inner_disp .^ 2, dims=2))) / h0
            if max_disp < dptol
                break
            end
        end
    end

    # Final triangulation
    t = triangulate_2d(triangulation, p)
    if size(t, 1) > 0
        pmid = (p[t[:, 1], :] .+ p[t[:, 2], :] .+ p[t[:, 3], :]) ./ 3
        d_mid = fd(pmid)
        t = t[vec(d_mid) .< -geps, :]
    end

    return p, t
end
