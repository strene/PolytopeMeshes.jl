"""
Conflict point removal utilities.
"""

"""
    remove_conflict_points(P1, P2, dist)

Remove points from P1 that are closer than the allowed distance to points in P2.

# Arguments
- `P1`: n×2 array of candidate points
- `P2`: m×2 array of reference points
- `dist`: m-element vector of minimum allowed distances from P2

# Returns
- `P1_filtered`: Points from P1 that are far enough from all P2 points
- `removed`: Boolean array indicating which P1 points were removed
"""
function remove_conflict_points(P1::Matrix{Float64}, P2::Matrix{Float64},
    dist::Vector{Float64})
    if isempty(P2) || size(P2, 1) == 0
        return copy(P1), falses(size(P1, 1))
    end

    n1 = size(P1, 1)
    n2 = size(P2, 1)
    removed = falses(n1)

    for i in 1:n1
        for j in 1:n2
            d = sqrt((P1[i, 1] - P2[j, 1])^2 + (P1[i, 2] - P2[j, 2])^2)
            if d < dist[j]
                removed[i] = true
                break
            end
        end
    end

    return P1[.!removed, :], removed
end
