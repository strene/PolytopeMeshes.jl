"""
Geometry utility functions for mesh generation.
"""

"""
    line_line_intersection(L1, L2)

Calculate the intersections of line segments L1 with line segments L2.

# Arguments
- `L1`: Matrix of line segments [n × 4], each row is [x1, y1, x2, y2]
- `L2`: Matrix of line segments [m × 4], each row is [x1, y1, x2, y2]

# Returns
- `X`: n × m matrix of x-coordinates of intersections
- `Y`: n × m matrix of y-coordinates of intersections
- `seg_int`: n × m boolean matrix, true where segments actually intersect
"""
function line_line_intersection(L1::Matrix{Float64}, L2::Matrix{Float64})
    TOL = 1000 * (maximum(L1) - minimum(L1)) * eps()
    n1 = size(L1, 1)
    n2 = size(L2, 1)

    X = zeros(n1, n2)
    Y = zeros(n1, n2)
    seg_int = falses(n1, n2)

    for i in 1:n1
        l1x1, l1y1, l1x2, l1y2 = L1[i, 1], L1[i, 2], L1[i, 3], L1[i, 4]
        A1 = l1y2 - l1y1
        B1 = l1x1 - l1x2
        C1 = A1 * l1x1 + B1 * l1y1

        for j in 1:n2
            l2x1, l2y1, l2x2, l2y2 = L2[j, 1], L2[j, 2], L2[j, 3], L2[j, 4]
            A2 = l2y2 - l2y1
            B2 = l2x1 - l2x2
            C2 = A2 * l2x1 + B2 * l2y1

            det_val = A1 * B2 - A2 * B1
            if abs(det_val) < eps()
                continue
            end

            x = (B2 * C1 - B1 * C2) / det_val
            y = (A1 * C2 - A2 * C1) / det_val

            X[i, j] = x
            Y[i, j] = y

            seg_int[i, j] = (min(l1x1, l1x2) - x <= TOL &&
                             x - max(l1x1, l1x2) <= TOL &&
                             min(l2x1, l2x2) - x <= TOL &&
                             x - max(l2x1, l2x2) <= TOL &&
                             min(l1y1, l1y2) - y <= TOL &&
                             y - max(l1y1, l1y2) <= TOL &&
                             min(l2y1, l2y2) - y <= TOL &&
                             y - max(l2y1, l2y2) <= TOL)
        end
    end

    return X, Y, seg_int
end

"""
    euclidean_distance(a, b)

Compute Euclidean distance between points (row vectors or matrices).
"""
function euclidean_distance(a::AbstractMatrix, b::AbstractMatrix)
    return sqrt.(sum((a .- b) .^ 2, dims=2))
end

function euclidean_distance(a::AbstractVector, b::AbstractVector)
    return norm(a - b)
end

"""
    point_in_polygon(px, py, poly)

Test if point (px, py) is inside polygon `poly` (n×2 matrix).
Uses the ray casting algorithm.
"""
function point_in_polygon(px::Float64, py::Float64, poly::Matrix{Float64})
    n = size(poly, 1)
    inside = false
    j = n
    for i in 1:n
        xi, yi = poly[i, 1], poly[i, 2]
        xj, yj = poly[j, 1], poly[j, 2]
        if ((yi > py) != (yj > py)) &&
           (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""
    points_in_polygon(pts, poly)

Test if each point in `pts` (n×2) is inside polygon `poly` (m×2).
"""
function points_in_polygon(pts::Matrix{Float64}, poly::Matrix{Float64})
    n = size(pts, 1)
    result = BitVector(undef, n)
    for i in 1:n
        result[i] = point_in_polygon(pts[i, 1], pts[i, 2], poly)
    end
    return result
end

"""
    insert_vec(A, B, id)

Insert rows of B into A at positions given by id.
Ported from MRST's insertVec.
"""
function insert_vec(A::AbstractMatrix, B::AbstractMatrix, id::AbstractVector{<:Integer})
    sorted_idx = sortperm(id)
    id_sorted = id[sorted_idx]
    B_sorted = B[sorted_idx, :]

    n_a = size(A, 1)
    n_b = size(B_sorted, 1)
    cols = size(A, 2)
    C = zeros(n_a + n_b, cols)

    # Positions where B elements go
    b_positions = id_sorted .+ collect(0:(n_b - 1))

    # Fill in B at designated positions
    a_idx = 1
    b_idx = 1
    for c_idx in 1:(n_a + n_b)
        if b_idx <= n_b && c_idx == b_positions[b_idx]
            C[c_idx, :] = B_sorted[b_idx, :]
            b_idx += 1
        else
            C[c_idx, :] = A[a_idx, :]
            a_idx += 1
        end
    end

    return C
end

function insert_vec(A::AbstractVector, B::AbstractVector, id::AbstractVector{<:Integer})
    sorted_idx = sortperm(id)
    id_sorted = id[sorted_idx]
    B_sorted = B[sorted_idx]

    n_a = length(A)
    n_b = length(B_sorted)
    C = zeros(eltype(A), n_a + n_b)

    b_positions = id_sorted .+ collect(0:(n_b - 1))

    a_idx = 1
    b_idx = 1
    for c_idx in 1:(n_a + n_b)
        if b_idx <= n_b && c_idx == b_positions[b_idx]
            C[c_idx] = B_sorted[b_idx]
            b_idx += 1
        else
            C[c_idx] = A[a_idx]
            a_idx += 1
        end
    end

    return C
end
