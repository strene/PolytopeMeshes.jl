# Default triangulation backend: the DistMesh algorithm (Per-Olof Persson and
# Gilbert Strang, "A Simple Mesh Generator in MATLAB", SIAM Review 46(2),
# 2004) from DistMesh.jl. The wrapper reproduces the behavior of the modified
# version in UPR (distmesh2d): the fixed points are kept in front of the
# returned point set, the iterations are capped, and the result is
# reproducible.

"""
    DistMeshBackend(; maxit = 200)

Triangulation backend using DistMesh.jl (the default). `maxit` caps the
number of force-equilibrium iterations, as in UPR.
"""
Base.@kwdef struct DistMeshBackend <: TriangulationBackend
    maxit::Int = 200
end

function triangulate_polygon(b::DistMeshBackend, polygon, size, fixed; rng = Random.default_rng(), min_size)
    fd = distance_function_polygon(polygon)
    lo, hi = bounding_box(polygon)
    return _distmesh(fd, x -> size(x) / min_size, min_size, (lo, hi), fixed; rng = rng, maxit = b.maxit)
end

# Smallest relative size requested inside the domain, sampled on the initial
# lattice of spacing h0.
function _min_size(fd, fh, h0, lo::Vec2, hi::Vec2)
    hmin = Inf
    geps = 0.001 * h0
    for x in lo[1]:h0:hi[1], y in lo[2]:(h0 * sqrt(3) / 2):hi[2]
        p = Vec2(x, y)
        if fd(p) < geps
            hmin = min(hmin, fh(p))
        end
    end
    return isfinite(hmin) && hmin > 0 ? hmin : 1.0
end

# Run DistMesh.jl. Returns the points (fixed points first, in the given
# order) and the triangles of the last retriangulation.
function _distmesh(fd, fh, h0, bbox, pfix; rng = Random.default_rng(), maxit = 200)
    lo, hi = Vec2(bbox[1]), Vec2(bbox[2])
    pfix = to_points(pfix, 2)
    pfix, _ = unique_points(pfix; tol = 1e-12 * norm(hi - lo))
    nfix = length(pfix)
    # DistMesh.jl normalizes the rejection probability by its maximum, while
    # UPR does not. Rescale so that the smallest requested size is h0, which
    # gives the same initial density as UPR.
    hmin = _min_size(fd, fh, h0, lo, hi)
    h = h0 * hmin
    hrel = x -> fh(Vec2(x)) / hmin
    dist = x -> fd(Vec2(x))
    # DistMesh.jl uses the task-local random number generator. Seed it in a
    # separate task to get reproducible grids without touching the caller's.
    seed = rand(rng, UInt64)
    msh = fetch(@async begin
        Random.seed!(seed)
        DistMesh.distmesh2d(dist, hrel, h, (Tuple(lo), Tuple(hi)), pfix; maxiter = maxit)
    end)
    p = Vec2.(msh.p)
    # DistMesh.jl prepends the fixed points, but may project those on the
    # boundary by round-off. Pin them to their exact positions.
    @assert length(p) >= nfix && all(i -> norm(p[i] - pfix[i]) <= 1e-8 * h, 1:nfix) "DistMesh.jl did not keep the fixed points first"
    p[1:nfix] .= pfix
    t = [(Int(T[1]), Int(T[2]), Int(T[3])) for T in msh.t]
    return p, t
end

"""
    distmesh_2d(fd, fh, h0, bbox, pfix; kwarg...)

Distribute points in the 2D domain given by the signed distance function
`fd(x)` (negative inside) with relative point spacing `fh(x)` and initial
spacing `h0`, using DistMesh.jl. `bbox = (lo, hi)` bounds the domain and
`pfix` are fixed points. The returned points start with `pfix` (in the same
order, with duplicates removed).

# Keyword arguments
- `rng = Random.default_rng()`: random number generator used to seed the
  initial rejection sampling.
- `maxit = 200`: maximum number of iterations (as in UPR).
"""
function distmesh_2d(fd, fh, h0, bbox, pfix = Vec2[]; rng = Random.default_rng(), maxit = 200)
    p, _ = _distmesh(fd, fh, h0, bbox, pfix; rng = rng, maxit = maxit)
    nfix = length(unique_points(to_points(pfix, 2); tol = 1e-12 * norm(Vec2(bbox[2]) - Vec2(bbox[1])))[1])
    # Remove free points outside the domain and duplicates
    keep = [i <= nfix || fd(p[i]) < 0.001 * h0 for i in eachindex(p)]
    up, _ = unique_points(p[keep]; tol = 1e-10 * norm(Vec2(bbox[2]) - Vec2(bbox[1])))
    return up
end

"""
    distance_function_rectangle(lo, hi)

Signed distance function of the rectangle `[lo, hi]` for [`distmesh_2d`](@ref).
"""
distance_function_rectangle(lo, hi) = x -> signed_distance_rectangle(x, lo, hi)

"""
    distance_function_polygon(poly)

Signed distance function of a polygon for [`distmesh_2d`](@ref).
"""
function distance_function_polygon(poly)
    v = to_points(poly, 2)
    return x -> signed_distance_polygon(x, v)
end
