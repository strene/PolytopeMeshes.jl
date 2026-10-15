# Polyline utilities: arc length parametrization, splitting at intersections
# (UPR: splitAtInt2D) and size-function driven point placement along lines
# (UPR: interLinePath, subdivideLineSegments, eqInterpret).

function arclength(line::AbstractVector{<:SVector})
    s = zeros(length(line))
    for i in 2:length(line)
        s[i] = s[i - 1] + norm(line[i] - line[i - 1])
    end
    return s
end

"""
    point_at(line, s, t)

Point at arc length `t` along `line` with cumulative arc lengths `s`.
"""
function point_at(line::AbstractVector{<:SVector}, s::AbstractVector, t::Real)
    n = length(line)
    t <= s[1] && return line[1]
    t >= s[end] && return line[end]
    k = searchsortedlast(s, t)
    k = clamp(k, 1, n - 1)
    ds = s[k + 1] - s[k]
    w = ds > 0 ? (t - s[k]) / ds : 0.0
    return line[k] + w * (line[k + 1] - line[k])
end

"""
    split_at_intersections(L1, L2 = [])

Split the 2D polylines `L1` at their mutual intersections and at their
intersections with the polylines `L2` (UPR: `splitAtInt2D`).

Returns a named tuple with fields
- `lines`: the split polylines,
- `self_cut`: for each piece, `2*(start is cut by L1) + (end is cut by L1)`,
- `other_cut`: the same flags for cuts by `L2`,
- `other_start`, `other_end`: index of the `L2` line cutting the piece at its
  start/end (`0` if none),
- `origin`: index into `L1` of the line each piece came from.

Lines with a single point are passed through unchanged.
"""
function split_at_intersections(L1, L2 = Vector{Vec2}[]; tol = 1e-10)
    L1 = to_polylines(L1, 2)
    L2 = to_polylines(L2, 2)
    lines = Vector{Vec2}[]
    self_cut = Int[]
    other_cut = Int[]
    other_start = Int[]
    other_end = Int[]
    origin = Int[]
    for (i, l) in enumerate(L1)
        if length(l) < 2
            push!(lines, copy(l)); push!(self_cut, 0); push!(other_cut, 0)
            push!(other_start, 0); push!(other_end, 0); push!(origin, i)
            continue
        end
        s = arclength(l)
        L = s[end]
        # Cuts as (arc length, kind, other line) with kind 1 = L1, 2 = L2
        cuts = Tuple{Float64, Int, Int}[]
        for (kind, others) in ((1, L1), (2, L2))
            for (j, m) in enumerate(others)
                (kind == 1 && j == i) && continue
                length(m) < 2 && continue
                for a in 1:(length(l) - 1), b in 1:(length(m) - 1)
                    r = segment_intersection(l[a], l[a + 1], m[b], m[b + 1]; tol = tol)
                    isnothing(r) && continue
                    t = s[a] + r[1] * (s[a + 1] - s[a])
                    push!(cuts, (t, kind, j))
                end
            end
        end
        sort!(cuts)
        # Merge cuts at (almost) the same location, preferring L1 cuts
        merged = Tuple{Float64, Int, Int}[]
        for c in cuts
            if !isempty(merged) && abs(c[1] - merged[end][1]) <= tol * max(L, 1.0)
                if c[2] < merged[end][2]
                    merged[end] = (merged[end][1], c[2], c[3])
                end
            else
                push!(merged, c)
            end
        end
        # Pieces between consecutive cut locations
        starts = [(0.0, 0, 0)]
        ends = Tuple{Float64, Int, Int}[]
        for c in merged
            if c[1] <= tol * L
                starts[1] = (0.0, c[2], c[3])
            elseif c[1] >= L * (1 - tol)
                # Cut at the end of the line, handled below
            else
                push!(ends, c)
                push!(starts, c)
            end
        end
        last_cut = (L, 0, 0)
        for c in merged
            if c[1] >= L * (1 - tol)
                last_cut = (L, c[2], c[3])
            end
        end
        push!(ends, last_cut)
        for (st, en) in zip(starts, ends)
            t0, t1 = st[1], en[1]
            t1 - t0 <= tol * L && continue
            piece = [point_at(l, s, t0)]
            for k in 2:(length(l) - 1)
                if t0 + tol * L < s[k] < t1 - tol * L
                    push!(piece, l[k])
                end
            end
            push!(piece, point_at(l, s, t1))
            push!(lines, piece)
            push!(self_cut, 2 * (st[2] == 1) + (en[2] == 1))
            push!(other_cut, 2 * (st[2] == 2) + (en[2] == 2))
            push!(other_start, st[2] == 2 ? st[3] : 0)
            push!(other_end, en[2] == 2 ? en[3] : 0)
            push!(origin, i)
        end
    end
    return (lines = lines, self_cut = self_cut, other_cut = other_cut,
        other_start = other_start, other_end = other_end, origin = origin)
end

# Evaluate a size function given either as a number or a function of a point
_eval_size(f::Real, x) = Float64(f)
_eval_size(f, x) = Float64(f(x))

"""
    interpolate_line_path(line, fh, ds; se = (0, 0), interpolate = false)

Place points along the polyline `line` so that the distance between
consecutive points follows the size function `fh(x)` (UPR: `interLinePath`).
The initial guess uses spacing `ds`. If `interpolate` is `false` the
vertices of the polyline are kept fixed; otherwise only the end points are.

`se = (a, b)`: if nonzero, the first (last) point is placed at distance
`a*fh` (`b*fh`) from the start (end) of the line and the end points
themselves are not returned.
"""
function interpolate_line_path(line::AbstractVector{SVector{D, Float64}}, fh, ds; se = (0.0, 0.0),
        interpolate = false, tol = 1e-4, maxit = 10000) where {D}
    n = length(line)
    n == 0 && return SVector{D, Float64}[]
    n == 1 && return copy(line)
    s_line = arclength(line)
    L = s_line[end]
    s0 = se[1] * ds
    e0 = L - se[2] * ds
    if e0 <= s0
        # Line too short for the requested offsets: a single point in the middle
        return [point_at(line, s_line, (s0 + e0) / 2)]
    end
    # Initial distribution in arc length with fixed flags
    if interpolate
        m = max(ceil(Int, (e0 - s0) / ds - 1e-10), 1)
        s = collect(range(s0, e0, length = m + 1))
        fixed = falses(length(s))
    else
        s = Float64[]
        fixed = Bool[]
        for k in 1:(n - 1)
            a = max(s_line[k], s0)
            b = min(s_line[k + 1], e0)
            b < a && continue
            m = max(ceil(Int, (b - a) / ds - 1e-10), 1)
            for (j, t) in enumerate(range(a, b, length = m + 1))
                if !isempty(s) && j == 1
                    continue
                end
                push!(s, t)
                push!(fixed, false)
            end
        end
        # Interior polyline vertices are fixed
        for k in 2:(n - 1)
            if s0 < s_line[k] < e0
                idx = argmin(abs.(s .- s_line[k]))
                s[idx] = s_line[k]
                fixed[idx] = true
            end
        end
    end
    # The end points of the line are always fixed
    if se[1] != 0
        pushfirst!(s, 0.0); pushfirst!(fixed, true)
    else
        s[1] = 0.0; fixed[1] = true
    end
    if se[2] != 0
        push!(s, L); push!(fixed, true)
    else
        s[end] = L; fixed[end] = true
    end
    count = 0
    while count < maxit
        count += 1
        np = length(s)
        d = diff(s)
        pmid = [point_at(line, s_line, (s[k] + s[k + 1]) / 2) for k in 1:(np - 1)]
        dw = [_eval_size(fh, p) for p in pmid]
        if se[1] != 0
            dw[1] *= se[1]
        end
        if se[2] != 0
            dw[end] *= se[2]
        end
        dd = d .- dw
        # Insert or remove points between consecutive fixed points
        F = findall(fixed)
        changes = Tuple{Int, Int}[]  # (+1 insert at bar, -1 remove node)
        for g in 1:(length(F) - 1)
            ix = F[g]:(F[g + 1] - 1)
            dsum = sum(@view dd[ix])
            if dsum > minimum(@view dw[ix])
                b = ix[argmax(@view dd[ix])]
                push!(changes, (1, b))
            elseif dsum < -maximum(@view dw[ix])
                b = ix[argmin(@view dd[ix])]
                node = b == first(ix) ? b + 1 : b
                if !fixed[node]
                    push!(changes, (-1, node))
                end
            end
        end
        if !isempty(changes)
            ins = Set(b for (k, b) in changes if k == 1)
            rem = Set(b for (k, b) in changes if k == -1)
            snew = Float64[]
            fnew = Bool[]
            for k in 1:np
                if !(k in rem)
                    push!(snew, s[k]); push!(fnew, fixed[k])
                end
                if k in ins
                    push!(snew, (s[k] + s[k + 1]) / 2); push!(fnew, false)
                end
            end
            s, fixed = snew, fnew
            continue
        end
        np <= 2 && break
        # Spring relaxation of the free points
        Fb = dw .- d
        Fn = Fb[1:(end - 1)] .- Fb[2:end]
        move = 0.2 .* Fn
        maxmove = 0.0
        for k in 2:(np - 1)
            if !fixed[k]
                # Limit the movement so that points never pass each other
                lim = 0.45 * min(d[k - 1], d[k])
                mv = clamp(move[k - 1], -lim, lim)
                s[k] += mv
                maxmove = max(maxmove, abs(mv))
            end
        end
        if maxmove < tol * ds
            break
        end
    end
    if count == maxit
        @warn "Point placement along line did not converge"
    end
    pts = [point_at(line, s_line, t) for t in s]
    if se[1] != 0
        popfirst!(pts)
    end
    if se[2] != 0
        pop!(pts)
    end
    return pts
end
