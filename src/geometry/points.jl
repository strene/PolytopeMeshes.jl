# Conversion of user input to vectors of static points and tolerance-based
# merging of points.

"""
    to_points(x, D)

Convert `x` to a `Vector{SVector{D, Float64}}`. Accepts an `n × D` matrix (MRST
convention), a vector of vectors/tuples, or a single point.
"""
function to_points(x::AbstractMatrix{<:Real}, D::Int)
    if size(x, 2) == D
        return [SVector{D, Float64}(x[i, :]) for i in axes(x, 1)]
    elseif size(x, 1) == D
        return [SVector{D, Float64}(x[:, i]) for i in axes(x, 2)]
    else
        throw(ArgumentError("Point matrix of size $(size(x)) is not compatible with dimension $D"))
    end
end
to_points(x::AbstractVector{<:Real}, D::Int) = [SVector{D, Float64}(x)]
to_points(x::NTuple{N, <:Real}, D::Int) where {N} = [SVector{D, Float64}(x)]
function to_points(x::AbstractVector, D::Int)
    out = Vector{SVector{D, Float64}}(undef, length(x))
    for (i, p) in enumerate(x)
        length(p) == D || throw(ArgumentError("Point $i has $(length(p)) coordinates, expected $D"))
        out[i] = SVector{D, Float64}(Tuple(p))
    end
    return out
end

"""
    to_polylines(x, D)

Convert a collection of polylines (each anything accepted by
[`to_points`](@ref)) to a vector of point vectors. A single polyline given as
a matrix is wrapped in a vector.
"""
to_polylines(x::AbstractMatrix{<:Real}, D::Int) = [to_points(x, D)]
to_polylines(x::Nothing, D::Int) = Vector{SVector{D, Float64}}[]
function to_polylines(x::AbstractVector, D::Int)
    isempty(x) && return Vector{SVector{D, Float64}}[]
    if first(x) isa Real
        return [to_points(x, D)]
    end
    if first(x) isa SVector{D} || (first(x) isa Union{Tuple, AbstractVector} && length(first(x)) == D && first(first(x)) isa Real)
        # A single polyline given as a vector of points
        return [to_points(x, D)]
    end
    return [to_points(l, D) for l in x]
end

function bounding_box(pts::AbstractVector{SVector{D, T}}) where {D, T}
    lo = reduce((a, b) -> min.(a, b), pts)
    hi = reduce((a, b) -> max.(a, b), pts)
    return lo, hi
end

"""
    PointMerger{D}(tol)

Spatial hash used to merge points closer than `tol`.
"""
struct PointMerger{D}
    tol::Float64
    points::Vector{SVector{D, Float64}}
    buckets::Dict{NTuple{D, Int}, Vector{Int}}
end
PointMerger{D}(tol) where {D} = PointMerger{D}(Float64(tol), SVector{D, Float64}[], Dict{NTuple{D, Int}, Vector{Int}}())

_bucket(pm::PointMerger{D}, x) where {D} = ntuple(i -> floor(Int, x[i] / pm.tol), D)

function _neighbor_offsets(::Val{D}) where {D}
    return vec(collect(Iterators.product(ntuple(_ -> -1:1, D)...)))
end

"""
    merge_point!(pm, x)

Return the index of a point within `pm.tol` of `x`, inserting `x` if none
exists.
"""
function merge_point!(pm::PointMerger{D}, x::SVector{D}) where {D}
    k = _bucket(pm, x)
    tol2 = pm.tol^2
    for off in _neighbor_offsets(Val(D))
        kk = k .+ off
        ids = get(pm.buckets, kk, nothing)
        isnothing(ids) && continue
        for id in ids
            if sum(abs2, pm.points[id] - x) <= tol2
                return id
            end
        end
    end
    push!(pm.points, x)
    id = length(pm.points)
    push!(get!(Vector{Int}, pm.buckets, k), id)
    return id
end

"""
    unique_points(pts; tol)

Merge points closer than `tol`, keeping the first occurrence. Returns
`(upts, map)` where `upts[map[i]]` is the representative of `pts[i]`.
"""
function unique_points(pts::AbstractVector{SVector{D, T}}; tol = 1e-10) where {D, T}
    pm = PointMerger{D}(tol)
    map = Vector{Int}(undef, length(pts))
    for (i, p) in enumerate(pts)
        map[i] = merge_point!(pm, SVector{D, Float64}(p))
    end
    return pm.points, map
end

"""
    min_distance(x, tree)

Distance from each point in `x` to the closest point in the KD-tree `tree`.
"""
function min_distance(x::AbstractVector{SVector{D, Float64}}, pts::AbstractVector{SVector{D, Float64}}) where {D}
    isempty(pts) && return fill(Inf, length(x))
    tree = KDTree(pts)
    _, d = nn(tree, x)
    return d
end
