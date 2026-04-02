"""
Split paths at intersections.
Ported from MRST's splitAtInt2D and splitLines.
"""

"""
    split_at_intersections_2d(L1, L2)

Split paths in L1 at all intersections with other paths in L1 and paths in L2.

# Arguments
- `L1`: Vector of n×2 matrices, each defining a piecewise linear path
- `L2`: Vector of n×2 matrices, each defining a piecewise linear path

# Returns
- `split_L1`: Vector of split path segments
- `L1_cut`: Cut indicators for L1-L1 intersections
- `L1_L2_cut`: Cut indicators for L1-L2 intersections
- `IC`: Index map from split paths back to original L1 paths
"""
function split_at_intersections_2d(
    L1::Vector{Matrix{Float64}},
    L2::Vector{Matrix{Float64}}
)
    if isempty(L1)
        return Matrix{Float64}[], Int[], Int[], Int[]
    end

    # First split L1 paths against each other
    split_results = Vector{Vector{Matrix{Float64}}}(undef, length(L1))
    tmp_cut = Vector{Vector{Int}}(undef, length(L1))
    IC = collect(1:length(L1))

    for i in 1:length(L1)
        others = vcat(L1[1:i-1], L1[i+1:end])
        splits, cuts = _split_lines([L1[i]], others)
        split_results[i] = splits
        tmp_cut[i] = cuts
    end

    # Flatten
    new_IC = Int[]
    all_splits = Matrix{Float64}[]
    all_tmp_cut = Int[]
    for i in 1:length(L1)
        for (s, c) in zip(split_results[i], tmp_cut[i])
            push!(all_splits, s)
            push!(all_tmp_cut, c)
            push!(new_IC, i)
        end
    end

    # Then split against L2
    if isempty(L2)
        return all_splits, all_tmp_cut, zeros(Int, length(all_splits)), new_IC
    end

    final_splits = Matrix{Float64}[]
    L1_cut = Int[]
    L1_L2_cut = Int[]
    final_IC = Int[]

    for i in 1:length(all_splits)
        splits, cuts = _split_lines([all_splits[i]], L2)
        n_new = length(splits)

        for (s, c) in zip(splits, cuts)
            push!(final_splits, s)
            push!(L1_L2_cut, c)
            push!(final_IC, new_IC[i])
        end

        # Propagate L1-L1 cut info
        tc = all_tmp_cut[i]
        if n_new == 1
            push!(L1_cut, tc)
        else
            for j in 1:n_new
                if j == 1
                    push!(L1_cut, tc >= 2 ? 2 : 0)
                elseif j == n_new
                    push!(L1_cut, (tc == 1 || tc == 3) ? 1 : 0)
                else
                    push!(L1_cut, 0)
                end
            end
        end
    end

    # Remove single-point lines
    keep = [size(s, 1) > 1 for s in final_splits]
    return final_splits[keep], L1_cut[keep], L1_L2_cut[keep], final_IC[keep]
end

"""
    _split_lines(L1, L2)

Split paths in L1 at intersections with paths in L2.
Returns split paths and cut indicators.
"""
function _split_lines(L1::Vector{Matrix{Float64}}, L2::Vector{Matrix{Float64}})
    split_lines = Matrix{Float64}[]
    is_cut = Int[]

    if isempty(L2)
        return copy(L1), zeros(Int, length(L1))
    end

    for i in 1:length(L1)
        l1 = L1[i]

        # Convert L2 to line segments [x1, y1, x2, y2]
        l2_segs = Vector{Vector{Float64}}()
        for l in L2
            for j in 1:(size(l, 1) - 1)
                push!(l2_segs, [l[j, 1], l[j, 2], l[j + 1, 1], l[j + 1, 2]])
            end
        end
        if isempty(l2_segs)
            push!(split_lines, l1)
            push!(is_cut, 0)
            continue
        end

        l2_mat = reduce(vcat, [s' for s in l2_segs])

        # Convert L1 to line segments
        l1_segs = zeros(size(l1, 1) - 1, 4)
        for j in 1:(size(l1, 1) - 1)
            l1_segs[j, :] = [l1[j, 1] l1[j, 2] l1[j + 1, 1] l1[j + 1, 2]]
        end

        X, Y, seg_int = line_line_intersection(l1_segs, l2_mat)

        # Find intersections
        int_points = Vector{Tuple{Float64,Float64,Int}}()  # (x, y, segment_index)
        for si in 1:size(seg_int, 1)
            for sj in 1:size(seg_int, 2)
                if seg_int[si, sj]
                    push!(int_points, (X[si, sj], Y[si, sj], si))
                end
            end
        end

        if isempty(int_points)
            push!(split_lines, l1)
            push!(is_cut, 0)
            continue
        end

        # Sort intersections by distance from start
        sort!(int_points, by=ip -> (ip[3], (ip[1] - l1[1, 1])^2 + (ip[2] - l1[1, 2])^2))

        # Remove duplicate intersection points
        unique_ints = [int_points[1]]
        for k in 2:length(int_points)
            prev = unique_ints[end]
            curr = int_points[k]
            if sqrt((curr[1] - prev[1])^2 + (curr[2] - prev[2])^2) > 1e-10
                push!(unique_ints, curr)
            end
        end

        # Split the line at intersection points
        _split_line_at_points!(split_lines, is_cut, l1, unique_ints)
    end

    return split_lines, is_cut
end

"""
Split a single line at given intersection points.
"""
function _split_line_at_points!(
    split_lines::Vector{Matrix{Float64}},
    is_cut::Vector{Int},
    line::Matrix{Float64},
    int_points::Vector{Tuple{Float64,Float64,Int}}
)
    current_path = [line[1, :]]
    start_cut = false
    n_splits = length(int_points)

    for (idx, (ix, iy, seg_idx)) in enumerate(int_points)
        int_pt = [ix, iy]

        # Add remaining original vertices up to this segment
        for j in 2:size(line, 1)
            if j - 1 <= seg_idx
                pt = line[j, :]
                # Check if we've passed the intersection
                if j - 1 == seg_idx
                    # Add the intersection point
                    if norm(current_path[end] - int_pt) > 1e-13
                        push!(current_path, int_pt)
                    end
                    # Save this segment
                    if length(current_path) >= 2
                        seg = reduce(vcat, [p' for p in current_path])
                        # Remove duplicate rows
                        seg = _unique_rows(seg)
                        if size(seg, 1) >= 2
                            push!(split_lines, seg)
                            cut_val = 0
                            if idx == 1 && start_cut
                                cut_val += 2
                            end
                            cut_val += 1  # end cut (intersection)
                            push!(is_cut, cut_val)
                        end
                    end
                    # Start new segment
                    current_path = [int_pt]
                    start_cut = true
                    break
                else
                    if isempty(current_path) || norm(current_path[end] - pt) > 1e-13
                        push!(current_path, pt)
                    end
                end
            end
        end
    end

    # Add remaining vertices after last intersection
    last_seg = int_points[end][3]
    for j in (last_seg + 1):size(line, 1)
        pt = line[j, :]
        if isempty(current_path) || norm(current_path[end] - pt) > 1e-13
            push!(current_path, pt)
        end
    end

    if length(current_path) >= 2
        seg = reduce(vcat, [p' for p in current_path])
        seg = _unique_rows(seg)
        if size(seg, 1) >= 2
            push!(split_lines, seg)
            push!(is_cut, start_cut ? 2 : 0)
        end
    end
end

function _unique_rows(A::Matrix{Float64}; tol::Float64 = 1e-13)
    if size(A, 1) <= 1
        return A
    end
    keep = [true]
    for i in 2:size(A, 1)
        if norm(A[i, :] - A[i - 1, :]) > tol
            push!(keep, true)
        else
            push!(keep, false)
        end
    end
    return A[keep, :]
end
