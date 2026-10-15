# Sites along cell constraints (wells, point constraints).

_per_line(x, n) = x isa AbstractVector ? (length(x) == n ? x : throw(ArgumentError("Expected $n values, got $(length(x))"))) : fill(x, n)

"""
    line_sites_2d(lines, ds; kwarg...)

Place sites along the 2D polylines `lines` so that the cells of the sites
follow the lines (UPR: `lineSites2D`). A line with a single point is a point
constraint.

# Keyword arguments
- `rho = ds`: absolute site spacing, a number or a function of position.
- `interpolate = false`: per line, see [`interpolate_line_path`](@ref).
- `se = (0, 0)`: per line offsets of the first/last site relative to `ds`.
- `cut = 0`: per line flag (`2*start + end`) telling if the line is cut by
  another cell constraint. Used to avoid protection sites at intersections.
- `protection_layer = false`: add sites on both sides of the lines.
- `protection_distance = ds/10`: per line distance (number or function of
  position) from the line to the protection sites.
- `line_index = 1:length(lines)`: index reported for the sites of each line.

Returns a named tuple with fields `sites`, `spacing`, `line_index`,
`protection_sites` and `protection_spacing`.
"""
function line_sites_2d(lines, ds; rho = ds, interpolate = false, se = (0.0, 0.0), cut = 0,
        protection_layer = false, protection_distance = ds / 10, line_index = nothing)
    lines = to_polylines(lines, 2)
    n = length(lines)
    interpolate = _per_line(interpolate, n)
    se = se isa Tuple ? fill(se, n) : _per_line(se, n)
    cut = _per_line(cut, n)
    protection_distance = _per_line(protection_distance, n)
    line_index = isnothing(line_index) ? collect(1:n) : line_index
    sites = Vec2[]
    spacing = Float64[]
    index = Int[]
    prot = Vec2[]
    prot_spacing = Float64[]
    for i in 1:n
        line = lines[i]
        if length(line) == 1
            p = copy(line)
            space = [Float64(_eval_size(rho, line[1]))]
        else
            p = interpolate_line_path(line, rho, ds; se = se[i], interpolate = interpolate[i])
            isempty(p) && continue
            if length(p) == 1
                space = [Float64(_eval_size(rho, p[1]))]
            else
                d = [norm(p[k + 1] - p[k]) for k in 1:(length(p) - 1)]
                space = [min(d[max(k - 1, 1)], d[min(k, length(d))]) for k in 1:length(p)]
            end
        end
        append!(sites, p)
        append!(spacing, space)
        append!(index, fill(line_index[i], length(p)))
        if protection_layer && length(line) > 1
            keep = collect(1:length(p))
            if cut[i] in (2, 3)
                keep = keep[2:end]
            end
            if cut[i] in (1, 3)
                keep = keep[1:(end - 1)]
            end
            isempty(keep) && continue
            pk = p[keep]
            if length(pk) == 1
                t = line[end] - pk[1]
                normals = [Vec2(t[2], -t[1]) / norm(t)]
            else
                normals = Vec2[]
                for k in 1:(length(pk) - 1)
                    t = pk[k + 1] - pk[k]
                    push!(normals, Vec2(t[2], -t[1]) / norm(t))
                end
                push!(normals, normals[end])
            end
            dist = [_eval_size(protection_distance[i], x) for x in pk]
            append!(prot, pk .+ normals .* dist)
            append!(prot, pk .- normals .* dist)
            append!(prot_spacing, dist)
            append!(prot_spacing, dist)
        end
    end
    # Remove duplicates (e.g. at intersections of lines)
    usites, map = unique_points(sites; tol = 1e-12 * max(ds, 1e-300))
    first_index = zeros(Int, length(usites))
    for (i, m) in enumerate(map)
        if first_index[m] == 0
            first_index[m] = i
        end
    end
    return (sites = usites, spacing = spacing[first_index], line_index = index[first_index],
        protection_sites = prot, protection_spacing = prot_spacing)
end

"""
    line_sites_3d(lines, rho)

Place sites along the 3D polylines `lines` with spacing given by `rho` (a
number, a function of position, or a vector with one per line) (UPR:
`lineSites3D`).

Returns a named tuple with fields `sites`, `spacing`, `line_index` and
`paths` (for each line, the consecutive sites as used by
[`line_sufficient_condition`](@ref)).
"""
function line_sites_3d(lines, rho)
    lines = to_polylines(lines, 3)
    n = length(lines)
    rho = rho isa AbstractVector ? _per_line(rho, n) : fill(rho, n)
    paths = Vector{Vec3}[]
    for i in 1:n
        l = lines[i]
        ds = _eval_size(rho[i], (l[1] + l[end]) / 2)
        p = interpolate_line_path(l, rho[i], ds)
        push!(paths, p)
    end
    allpts = reduce(vcat, paths; init = Vec3[])
    index = reduce(vcat, [fill(i, length(paths[i])) for i in 1:n]; init = Int[])
    L = 1.0
    if !isempty(allpts)
        lo, hi = bounding_box(allpts)
        L = max(norm(hi - lo), 1e-300)
    end
    usites, map = unique_points(allpts; tol = 1e-13 * L)
    first_index = zeros(Int, length(usites))
    for (i, m) in enumerate(map)
        if first_index[m] == 0
            first_index[m] = i
        end
    end
    spacing = fill(Inf, length(usites))
    offset = 0
    for p in paths
        for k in eachindex(p)
            d = Inf
            k > 1 && (d = min(d, norm(p[k] - p[k - 1])))
            k < length(p) && (d = min(d, norm(p[k + 1] - p[k])))
            id = map[offset + k]
            spacing[id] = min(spacing[id], d)
        end
        offset += length(p)
    end
    spacing[.!isfinite.(spacing)] .= 0.0
    return (sites = usites, spacing = spacing, line_index = index[first_index], paths = paths)
end
