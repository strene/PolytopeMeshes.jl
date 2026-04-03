"""
Construct a 2D clipped Voronoi mesh.
Ported from MRST's clippedPebi2D.
"""

"""
    clipped_voronoi_2d(p, bnd; triangulation=default_triangulation_backend())

Construct a 2D Voronoi mesh by computing the Voronoi diagram of points `p`
clipped to the polygon `bnd`.

# Arguments
- `p`: n×2 array of Voronoi site coordinates
- `bnd`: k×2 array of polygon boundary vertices (clockwise or counter-clockwise)

# Keyword Arguments
- `triangulation`: A [`TriangulationBackend`](@ref) instance to use for Delaunay
  triangulation (default: [`default_triangulation_backend()`](@ref))

# Returns
- `G::UnstructuredMesh`: Valid mesh definition

# Example
```julia
p = rand(30, 2)
bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
G = clipped_voronoi_2d(p, bnd)
```
"""
function clipped_voronoi_2d(p::Matrix{Float64}, bnd::Matrix{Float64};
    triangulation::TriangulationBackend = default_triangulation_backend())
    n_sites = size(p, 1)

    # Compute Delaunay triangulation
    tri = triangulate_2d(triangulation, p)
    edges = _get_edges(tri)

    # For each Voronoi cell, clip against boundary
    all_vertices = Vector{Float64}[]
    cell_vertex_indices = Vector{Vector{Int}}(undef, n_sites)

    for s in 1:n_sites
        # Find all edges connected to site s
        neighbor_edges = findall(i -> edges[i, 1] == s || edges[i, 2] == s, 1:size(edges, 1))

        # Get neighbor sites
        neighbors = Int[]
        for ei in neighbor_edges
            if edges[ei, 1] == s
                push!(neighbors, edges[ei, 2])
            else
                push!(neighbors, edges[ei, 1])
            end
        end

        if isempty(neighbors)
            cell_vertex_indices[s] = Int[]
            continue
        end

        # Compute bisector normals and midpoints
        n_neigh = length(neighbors)
        normals = zeros(n_neigh, 2)
        midpoints = zeros(n_neigh, 2)

        for j in 1:n_neigh
            nb = neighbors[j]
            diff_vec = p[nb, :] .- p[s, :]
            len = norm(diff_vec)
            if len > 0
                normals[j, :] = diff_vec ./ len
            end
            midpoints[j, :] = (p[nb, :] .+ p[s, :]) ./ 2
        end

        # Clip boundary polygon against all bisector half-planes
        cell_poly = clip_polygon(bnd, normals, midpoints)

        if isempty(cell_poly) || size(cell_poly, 1) < 3
            cell_vertex_indices[s] = Int[]
            continue
        end

        # Round to avoid floating point duplicates
        cell_poly = round.(cell_poly .* 1e6) ./ 1e6

        # Remove duplicate vertices (maintaining order)
        unique_poly = _unique_rows_stable(cell_poly)

        if size(unique_poly, 1) < 3
            cell_vertex_indices[s] = Int[]
            continue
        end

        # Store vertices
        indices = Int[]
        for j in 1:size(unique_poly, 1)
            v = unique_poly[j, :]
            # Check if this vertex already exists
            found = false
            for k in 1:length(all_vertices)
                if norm(all_vertices[k] - v) < 1e-6
                    push!(indices, k)
                    found = true
                    break
                end
            end
            if !found
                push!(all_vertices, v)
                push!(indices, length(all_vertices))
            end
        end

        cell_vertex_indices[s] = indices
    end

    # Build mesh structure
    if isempty(all_vertices)
        error("No valid cells generated")
    end

    V = reduce(vcat, [v' for v in all_vertices])
    n_nodes = size(V, 1)

    # Build faces from cell edges
    face_map = Dict{Tuple{Int,Int},Int}()
    face_nodes_list = Vector{Vector{Int}}()
    cell_faces_list = Vector{Vector{Int}}()
    face_neighbors = Vector{Vector{Int}}()

    for s in 1:n_sites
        idx = cell_vertex_indices[s]
        if isempty(idx)
            push!(cell_faces_list, Int[])
            continue
        end

        nv = length(idx)
        cf = Int[]

        for j in 1:nv
            k = mod(j, nv) + 1
            n1, n2 = idx[j], idx[k]
            edge_key = n1 < n2 ? (n1, n2) : (n2, n1)

            if haskey(face_map, edge_key)
                fi = face_map[edge_key]
                push!(cf, fi)
                # Add second neighbor
                face_neighbors[fi][2] = s
            else
                push!(face_nodes_list, [n1, n2])
                fi = length(face_nodes_list)
                face_map[edge_key] = fi
                push!(cf, fi)
                push!(face_neighbors, [s, 0])
            end
        end

        push!(cell_faces_list, cf)
    end

    # Build arrays
    n_cells = n_sites
    n_faces = length(face_nodes_list)

    # cells.faces and cells.facePos
    cell_faces = Int[]
    cell_facePos = [1]
    for cf in cell_faces_list
        append!(cell_faces, cf)
        push!(cell_facePos, cell_facePos[end] + length(cf))
    end

    # faces.nodes and faces.nodePos
    face_nodes = Int[]
    face_nodePos = [1]
    for fn in face_nodes_list
        append!(face_nodes, fn)
        push!(face_nodePos, face_nodePos[end] + length(fn))
    end

    # faces.neighbors
    fn_mat = zeros(Int, n_faces, 2)
    for i in 1:n_faces
        fn_mat[i, 1] = face_neighbors[i][1]
        fn_mat[i, 2] = face_neighbors[i][2]
    end

    G = UnstructuredMesh(V, cell_faces, cell_facePos, face_nodes, face_nodePos,
        fn_mat; meshdim=2)

    # Sort edges for consistent orientation
    sort_edges!(G)

    return G
end

"""
    _delaunay_2d(points)

Compute 2D Delaunay triangulation using the Bowyer-Watson algorithm.
Returns n×3 matrix of triangle vertex indices.
"""
function _delaunay_2d(points::Matrix{Float64})
    n = size(points, 1)
    if n < 3
        return zeros(Int, 0, 3)
    end

    # Create super triangle well outside the point set
    min_x = minimum(points[:, 1])
    max_x = maximum(points[:, 1])
    min_y = minimum(points[:, 2])
    max_y = maximum(points[:, 2])

    dx = max_x - min_x + 1e-10
    dy = max_y - min_y + 1e-10
    d_max = max(dx, dy)
    mid_x = (min_x + max_x) / 2
    mid_y = (min_y + max_y) / 2

    # Super triangle vertices
    p1 = [mid_x - 20 * d_max, mid_y - d_max]
    p2 = [mid_x, mid_y + 20 * d_max]
    p3 = [mid_x + 20 * d_max, mid_y - d_max]

    aug_pts = vcat(points, p1', p2', p3')
    st1, st2, st3 = n + 1, n + 2, n + 3

    # Use Vector for triangles (more predictable than Set for iteration)
    triangles = Vector{Tuple{Int,Int,Int}}()
    push!(triangles, _orient_tri(aug_pts, st1, st2, st3))

    for i in 1:n
        pt = vec(aug_pts[i, :])

        # Find bad triangles
        bad_idx = Int[]
        for (ti, tri) in enumerate(triangles)
            if _in_circumcircle(aug_pts, tri, pt)
                push!(bad_idx, ti)
            end
        end

        # Find boundary polygon (edges not shared by two bad triangles)
        boundary = Vector{Tuple{Int,Int}}()
        for ti in bad_idx
            tri = triangles[ti]
            edges_tri = [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
            for edge in edges_tri
                shared = false
                for tj in bad_idx
                    if tj == ti
                        continue
                    end
                    otri = triangles[tj]
                    oedges = [(otri[1], otri[2]), (otri[2], otri[3]), (otri[3], otri[1])]
                    for oe in oedges
                        if (edge[1] == oe[2] && edge[2] == oe[1])
                            shared = true
                            break
                        end
                    end
                    if shared
                        break
                    end
                end
                if !shared
                    push!(boundary, edge)
                end
            end
        end

        # Remove bad triangles (in reverse order to keep indices valid)
        sort!(bad_idx, rev=true)
        for ti in bad_idx
            deleteat!(triangles, ti)
        end

        # Create new triangles from boundary edges to new point
        for edge in boundary
            new_tri = _orient_tri(aug_pts, i, edge[1], edge[2])
            push!(triangles, new_tri)
        end
    end

    # Remove triangles that use super triangle vertices
    result = Vector{Vector{Int}}()
    for tri in triangles
        if tri[1] <= n && tri[2] <= n && tri[3] <= n
            push!(result, [tri[1], tri[2], tri[3]])
        end
    end

    if isempty(result)
        return zeros(Int, 0, 3)
    end

    return reduce(vcat, [r' for r in result])
end

# BowyerWatson backend method — defined here (not in triangulation.jl) because
# it depends on _delaunay_2d which is defined above in this file.
function triangulate_2d(::BowyerWatson, points::Matrix{Float64})
    return _delaunay_2d(points)
end

"""
Orient triangle vertices counter-clockwise.
"""
function _orient_tri(pts::Matrix{Float64}, i::Int, j::Int, k::Int)
    ax, ay = pts[i, 1], pts[i, 2]
    bx, by = pts[j, 1], pts[j, 2]
    cx, cy = pts[k, 1], pts[k, 2]
    cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if cross >= 0
        return (i, j, k)
    else
        return (i, k, j)
    end
end

"""
Test if point is inside the circumcircle of a triangle.
Uses the standard determinant-based incircle test.
Triangle vertices must be counter-clockwise oriented.
"""
function _in_circumcircle(pts::Matrix{Float64}, tri::Tuple{Int,Int,Int}, pt::Vector{Float64})
    ax = pts[tri[1], 1] - pt[1]
    ay = pts[tri[1], 2] - pt[2]
    bx = pts[tri[2], 1] - pt[1]
    by = pts[tri[2], 2] - pt[2]
    cx = pts[tri[3], 1] - pt[1]
    cy = pts[tri[3], 2] - pt[2]

    det = ax * (by * (cx^2 + cy^2) - cy * (bx^2 + by^2)) -
          ay * (bx * (cx^2 + cy^2) - cx * (bx^2 + by^2)) +
          (ax^2 + ay^2) * (bx * cy - by * cx)

    return det > 0
end

"""
Get unique edges from a triangulation.
"""
function _get_edges(tri::Matrix{Int})
    if size(tri, 1) == 0
        return zeros(Int, 0, 2)
    end

    edge_set = Set{Tuple{Int,Int}}()
    for i in 1:size(tri, 1)
        for (a, b) in [(tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 3]), (tri[i, 3], tri[i, 1])]
            e = a < b ? (a, b) : (b, a)
            push!(edge_set, e)
        end
    end

    edges = zeros(Int, length(edge_set), 2)
    for (i, (a, b)) in enumerate(edge_set)
        edges[i, :] = [a, b]
    end

    return edges
end

"""
Remove duplicate rows while maintaining order.
"""
function _unique_rows_stable(A::Matrix{Float64})
    seen = Set{Tuple{Float64,Float64}}()
    keep = Int[]
    for i in 1:size(A, 1)
        key = (A[i, 1], A[i, 2])
        if !(key in seen)
            push!(seen, key)
            push!(keep, i)
        end
    end
    return A[keep, :]
end
