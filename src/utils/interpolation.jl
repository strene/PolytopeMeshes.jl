"""
Line path interpolation utilities.
Ported from MRST's interLinePath, eqInterpret, subdivideLineSegments.
"""

"""
    eq_interpret(path, dt, se_ptn)

Interpolate a path with equidistant points.

# Arguments
- `path`: n×2 array of path coordinates
- `dt`: Distance between interpolation points
- `se_ptn`: [start, end] trim pattern (fraction of dt)

# Returns
- `new_points`: Interpolated points
- `dt`: Actual distance between points
"""
function eq_interpret(path::Matrix{Float64}, dt::Float64, se_ptn::Vector{Float64})
    lines_dist = sqrt.(sum(diff(path, dims=1) .^ 2, dims=2))
    lines_dist = [0.0; vec(lines_dist)]
    cum_dist = cumsum(lines_dist)

    s = dt * se_ptn[1]
    e = dt * se_ptn[2]
    total_len = cum_dist[end] - s - e
    if total_len <= 0
        return zeros(0, 2), dt
    end
    dt_new = total_len / ceil(total_len / dt)

    new_points_loc = collect(s:dt_new:(cum_dist[end] - e))
    # Make sure the last point is included
    if length(new_points_loc) == 0
        return zeros(0, 2), dt_new
    end

    new_points = _interp1(cum_dist, path, new_points_loc)
    return new_points, dt_new
end

"""
    subdivide_line_segments(path, dt, se_ptn)

Subdivide line segments with almost equidistant points, respecting input vertices.

# Arguments
- `path`: n×2 array of path coordinates
- `dt`: Default distance for subdivision
- `se_ptn`: Start/end trim pattern

# Returns
- `pts`: New points
- `orig_ix`: Indices of original path vertices among new points
- `pts_used`: Indices of original path vertices actually used
"""
function subdivide_line_segments(path::Matrix{Float64}, dt::Float64, se_ptn::Vector{Float64})
    s_dist = vec(sqrt.(sum(diff(path, dims=1) .^ 2, dims=2)))
    c_coor = [0.0; cumsum(s_dist)]
    split_counts = ceil.(Int, s_dist ./ dt)
    n_total = sum(split_counts) + 1
    n_coor = fill(NaN, n_total)

    s_pt = se_ptn[1] * dt
    e_pt = c_coor[end] - se_ptn[2] * dt
    n = 1
    flag = falses(length(s_dist) + 1)

    for i in 1:length(s_dist)
        if c_coor[i + 1] < s_pt || c_coor[i] > e_pt
            split_counts[i] = 0
            continue
        end
        ds = min(c_coor[i + 1], e_pt) - max(c_coor[i], s_pt)
        split_counts[i] = ceil(Int, ds / dt)
        start_val = max(c_coor[i], s_pt)
        step = ds / split_counts[i]
        for k in 0:split_counts[i]
            n_coor[n + k] = start_val + k * step
        end
        n += split_counts[i]
        flag[i] = true
        flag[i + 1] = true
    end

    n_coor = filter(!isnan, n_coor)
    valid_splits = filter(x -> x > 0, split_counts)

    orig_ix = [1; cumsum(valid_splits) .+ 1]
    if se_ptn[1] > 0
        orig_ix = orig_ix[2:end]
    end
    if se_ptn[2] > 0
        orig_ix = orig_ix[1:end-1]
    end

    pts = _interp1(c_coor, path, n_coor)
    pts_used = findall(flag)
    return pts, orig_ix, pts_used
end

"""
    inter_line_path(line, fh, line_dist, se_ptn, interpol)

Interpolate a line path with adaptive spacing.

# Arguments
- `line`: n×2 coordinates of the line path (ordered)
- `fh`: Distance function handle `fh(midpoints) -> distances`
- `line_dist`: Base distance between interpolation points
- `se_ptn`: Start/end trim pattern
- `interpol`: If true, use equidistant interpolation

# Returns
- `p`: Interpolated points (n×2 matrix)
"""
function inter_line_path(
    line::Matrix{Float64},
    fh::Function,
    line_dist::Float64,
    se_ptn::Vector{Float64},
    interpol::Bool
)
    TOL = 1e-4
    max_it = 10000

    if size(line, 1) == 1
        return copy(line)
    end
    if size(line, 1) == 0
        return zeros(0, 2)
    end

    if interpol
        p, _ = eq_interpret(line, line_dist, se_ptn)
        flag = falses(size(p, 1))
        pts_used = [1, size(line, 1)]
    else
        p, orig_ix, pts_used = subdivide_line_segments(line, line_dist, se_ptn)
        flag = falses(size(p, 1))
        for ix in orig_ix
            if ix >= 1 && ix <= length(flag)
                flag[ix] = true
            end
        end
        pts_used[[1, end]] = [1, size(line, 1)]
    end

    if size(p, 1) == 0
        return zeros(0, 2)
    end

    # Add auxiliary points
    if se_ptn[1] != 0
        p = vcat(line[1:1, :], p)
        flag = vcat([true], flag)
    end
    if se_ptn[2] != 0
        p = vcat(p, line[end:end, :])
        flag = vcat(flag, [true])
    end

    count = 0
    while count < max_it
        count += 1

        # Calculate distances and wanted distances
        d = _dist_along_line(line, p)
        pmid = (p[1:end-1, :] .+ p[2:end, :]) ./ 2
        dw = fh(pmid)

        if se_ptn[1] != 0
            dw[1] *= se_ptn[1]
        end
        if se_ptn[2] != 0
            dw[end] *= se_ptn[2]
        end

        # Determine whether to insert or remove points
        dd = d .- dw
        fix_indices = [0; findall(flag[2:end-1]); length(d)]
        insert_arr = zeros(Int, length(fix_indices) - 1)

        for i in 1:length(insert_arr)
            ix = (fix_indices[i] + 1):fix_indices[i + 1]
            if isempty(ix)
                continue
            end
            dsum = sum(dd[ix])
            if dsum > minimum(dw[ix])
                _, id = findmax(dd[ix])
                insert_arr[i] = id + fix_indices[i]
            elseif dsum < -maximum(dw[ix])
                _, id = findmin(dd[ix])
                id = id + fix_indices[i]
                if id == 1
                    id += 1
                end
                if !flag[id]
                    insert_arr[i] = -id
                end
            end
        end

        # Insert or remove points
        if any(insert_arr .!= 0)
            offset = 0
            for i in 1:length(insert_arr)
                if insert_arr[i] > 0
                    id = insert_arr[i] + offset
                    p = vcat(p[1:id, :], pmid[insert_arr[i]:insert_arr[i], :], p[id+1:end, :])
                    flag = vcat(flag[1:id], [false], flag[id+1:end])
                    offset += 1
                elseif insert_arr[i] < 0
                    id = -insert_arr[i] + offset
                    p = vcat(p[1:id-1, :], p[id+1:end, :])
                    flag = vcat(flag[1:id-1], flag[id+1:end])
                    offset -= 1
                end
            end
            continue
        end

        # If only external nodes
        if size(p, 1) <= 2
            break
        end

        # Move points based on desired length
        Fb = dw .- d
        Fn = Fb[1:end-1] .- Fb[2:end]
        move_node = Fn .* 0.2
        d .+= vcat([move_node[1]], move_node[2:end] .- move_node[1:end-1], [-move_node[end]])

        p = _interp_line(line, d, findall(flag), pts_used)

        # Terminate if nodes have moved less than TOL
        movable = .!flag[2:end-1]
        if all(abs.(move_node[movable]) .< TOL * line_dist)
            break
        end
    end

    # Remove auxiliary points
    if se_ptn[1] != 0
        p = p[2:end, :]
    end
    if se_ptn[2] != 0
        p = p[1:end-1, :]
    end

    return p
end

"""
    _dist_along_line(line, p)

Calculate distances between consecutive interpolation points along the line.
"""
function _dist_along_line(line::Matrix{Float64}, p::Matrix{Float64})
    TOL = 50 * eps()
    N = size(p, 1)
    d = zeros(N - 1)
    joint_dist = 0.0

    for i in 1:(size(line, 1) - 1)
        # Find which points of p lie on this line segment
        seg_start = line[i, :]
        seg_end = line[i + 1, :]
        seg_len = norm(seg_end - seg_start)

        indx = Int[]
        for k in 1:size(p, 1)
            dist_a = norm(p[k, :] - seg_start) + norm(seg_end - p[k, :])
            dist_b = seg_len
            if abs(dist_a - dist_b) < TOL * dist_b
                push!(indx, k)
            end
        end

        if isempty(indx)
            joint_dist += seg_len
            continue
        end

        if length(indx) >= 2
            for k in 1:(length(indx) - 1)
                d[indx[k]] = norm(p[indx[k], :] - p[indx[k + 1], :])
            end
        end

        if indx[1] > 1 && norm(seg_start - p[indx[1], :]) > TOL
            d[indx[1] - 1] = joint_dist + norm(seg_start - p[indx[1], :])
        end
        joint_dist = norm(p[indx[end], :] - seg_end)
    end

    return d
end

"""
    _interp_line(path, dt, fix, pts_used)

Interpolate along a path given distances between consecutive points.
"""
function _interp_line(path::Matrix{Float64}, dt::Vector{Float64},
    fix::Vector{Int}, pts_used::Vector{Int})
    dist_s = vec(sqrt.(sum(diff(path, dims=1) .^ 2, dims=2)))
    t = [0.0; cumsum(dist_s)]

    new_pts_eval = [0.0; cumsum(dt)]
    if length(fix) > 2
        for (i, f) in enumerate(fix)
            if i <= length(pts_used)
                new_pts_eval[f] = t[pts_used[i]]
            end
        end
    else
        new_pts_eval[end] = t[end]
    end

    return _interp1(t, path, new_pts_eval)
end

"""
    _interp1(x, y, xi)

1D interpolation (linear). Equivalent to MATLAB's interp1.
"""
function _interp1(x::Vector{Float64}, y::Matrix{Float64}, xi::Vector{Float64})
    n = length(x)
    m = size(y, 2)
    ni = length(xi)
    yi = zeros(ni, m)

    for i in 1:ni
        xv = clamp(xi[i], x[1], x[end])
        # Find interval
        idx = searchsortedlast(x, xv)
        if idx < 1
            idx = 1
        end
        if idx >= n
            idx = n - 1
        end
        # Linear interpolation
        t = (xv - x[idx]) / (x[idx + 1] - x[idx])
        if isnan(t) || isinf(t)
            t = 0.0
        end
        yi[i, :] = y[idx, :] .+ t .* (y[idx + 1, :] .- y[idx, :])
    end

    return yi
end
