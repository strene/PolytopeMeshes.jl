# Centroidal Voronoi (CPG) grids by Lloyd iterations. UPR (CPG2D/CPG3D)
# minimizes the CPG energy with L-BFGS; the gradient of that energy with
# respect to a site is 2*m*(x - centroid), and Lloyd's method is the
# corresponding fixed-point iteration.

# Degree-2 quadrature on the sub-simplices of the cells for the density
# weighted mass and centroid.
function _density_centroids(m::PolytopalMesh{D}, rho) where {D}
    geo = mesh_geometry(m)
    nc = number_of_cells(m)
    T = SVector{D, Float64}
    mass = zeros(nc)
    first_moment = zeros(T, nc)
    a, b = 0.5854101966249685, 0.1381966011250105
    for c in 1:nc
        x0 = geo.cell_centroids[c]
        for f in cell_faces(m, c)
            sgn = m.face_neighbors[f][1] == c ? 1.0 : -1.0
            nodes = face_nodes(m, f)
            if D == 2
                p1, p2 = m.node_points[nodes[1]], m.node_points[nodes[2]]
                vol = sgn * cross2(p1 - x0, p2 - x0) / 2
                for q in ((x0 + p1) / 2, (p1 + p2) / 2, (p2 + x0) / 2)
                    w = vol / 3 * rho(q)
                    mass[c] += w
                    first_moment[c] += w * q
                end
            else
                xf = geo.face_centroids[f]
                k = length(nodes)
                for i in 1:k
                    p1 = m.node_points[nodes[i]]
                    p2 = m.node_points[nodes[mod1(i + 1, k)]]
                    vol = sgn * dot(cross(p1 - xf, p2 - xf), xf - x0) / 6
                    v = (x0, xf, p1, p2)
                    for j in 1:4
                        q = b * (v[1] + v[2] + v[3] + v[4]) + (a - b) * v[j]
                        w = vol / 4 * rho(q)
                        mass[c] += w
                        first_moment[c] += w * q
                    end
                end
            end
        end
    end
    return [mass[c] > 0 ? first_moment[c] / mass[c] : geo.cell_centroids[c] for c in 1:nc]
end

"""
    centroidal_voronoi_mesh(sites, domain; kwarg...)

Optimize the positions of `sites` towards a centroidal Voronoi diagram in
`domain` by Lloyd iterations, and return `(mesh, sites)` (UPR: `CPG2D`,
`CPG3D`).

# Keyword arguments
- `fixed = 0`: number of leading sites that are kept fixed, or a boolean mask.
- `density = nothing`: density function `rho(x)`; denser regions get
  smaller cells.
- `maxit = 40`: maximum number of iterations.
- `tol = 1e-6`: stop when no site moves more than `tol` times the domain
  size.

Remaining keyword arguments are passed on to [`voronoi_mesh`](@ref) for the
final mesh.
"""
function centroidal_voronoi_mesh(sites, domain::AbstractDomain{D}; fixed = 0, density = nothing, maxit = 40, tol = 1e-6, kwarg...) where {D}
    x = to_points(sites, D)
    isfixed = fixed isa Integer ? [i <= fixed for i in eachindex(x)] : collect(Bool, fixed)
    L = domain_diameter(domain)
    for it in 1:maxit
        m, map = voronoi_mesh(x, domain; return_site_map = true)
        cents = isnothing(density) ? mesh_geometry(m).cell_centroids : _density_centroids(m, density)
        maxmove = 0.0
        for (c, i) in enumerate(map)
            isfixed[i] && continue
            maxmove = max(maxmove, norm(cents[c] - x[i]))
            x[i] = cents[c]
        end
        maxmove < tol * L && break
    end
    return voronoi_mesh(x, domain; kwarg...), x
end
