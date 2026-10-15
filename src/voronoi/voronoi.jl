# Restricted Voronoi diagrams computed cell by cell: the cell of a site is the
# domain clipped by the bisectors of the site and its nearest neighbors. The
# clipping stops when the next neighbor is further away than twice the
# distance to the farthest cell vertex ("security radius").

"""
    voronoi_mesh(sites, domain; kwarg...)

Compute the Voronoi diagram of `sites` restricted to `domain` and return it as
a [`PolytopalMesh`](@ref). Sites outside the domain are ignored and duplicate
sites (closer than `merge_tol`) are merged.

# Keyword arguments
- `cell_tags`: tag for each site, transferred to the cells.
- `face_tag`: function `(i, j) -> Int` giving the tag of the face between the
  cells of sites `i` and `j` (indices into `sites`).
- `tol`: tolerance used when clipping (default `1e-10` times the domain size).
- `merge_tol`: tolerance for merging vertices (default `1e-8` times the domain
  size).
- `return_site_map`: also return the index of the site of each cell.
"""
function voronoi_mesh(sites, domain::AbstractDomain{D};
        cell_tags = nothing,
        face_tag = nothing,
        tol = nothing,
        merge_tol = nothing,
        return_site_map = false
    ) where {D}
    pts = to_points(sites, D)
    L = domain_diameter(domain)
    tol = Float64(something(tol, 1e-10 * L))
    merge_tol = Float64(something(merge_tol, 1e-8 * L))
    # Keep unique sites inside the domain
    keep = Int[]
    pm = PointMerger{D}(merge_tol)
    ndup = 0
    for (i, p) in enumerate(pts)
        inside(domain, p; tol = tol) || continue
        n0 = length(pm.points)
        id = merge_point!(pm, p)
        if id > n0
            push!(keep, i)
        else
            ndup += 1
        end
    end
    if ndup > 0
        @warn "Removed $ndup duplicate sites"
    end
    x = pts[keep]
    cells = _voronoi_cells(x, domain, tol)
    mesh_data = _assemble_cells(cells, merge_tol, Val(D))
    nc = length(x)
    ctags = isnothing(cell_tags) ? zeros(Int, nc) : Int[cell_tags[keep[c]] for c in 1:nc]
    nf = length(mesh_data.neighbors)
    ftags = zeros(Int, nf)
    if !isnothing(face_tag)
        for f in 1:nf
            l, r = mesh_data.neighbors[f]
            if r > 0
                ftags[f] = face_tag(keep[l], keep[r])
            end
        end
    end
    mesh = PolytopalMesh{D}(
        mesh_data.nodes,
        mesh_data.faces_nodes,
        mesh_data.faces_nodespos,
        mesh_data.neighbors,
        mesh_data.cells_faces,
        mesh_data.cells_facepos,
        x,
        ctags,
        ftags,
        mesh_data.boundary_tags
    )
    if return_site_map
        return mesh, keep
    end
    return mesh
end

function _voronoi_cells(x::Vector{SVector{D, Float64}}, domain, tol) where {D}
    n = length(x)
    P0 = initial_polytope(domain)
    cells = Vector{typeof(P0)}(undef, n)
    n == 0 && return cells
    if n == 1
        cells[1] = P0
        return cells
    end
    tree = KDTree(x)
    k0 = min(n, D == 2 ? 16 : 40)
    Threads.@threads for i in 1:n
        cells[i] = _clip_cell(i, x, tree, P0, k0, tol)
    end
    return cells
end

function _clip_cell(i, x::Vector{SVector{D, Float64}}, tree, P0, k0, tol) where {D}
    n = length(x)
    xi = x[i]
    P = P0
    R = max_radius(P, xi)
    k = k0
    processed = Set{Int}()
    while true
        idxs, dists = knn(tree, xi, k, true)
        finished = false
        for m in eachindex(idxs)
            j = idxs[m]
            j == i && continue
            if dists[m] > 2 * R + tol
                finished = true
                break
            end
            j in processed && continue
            push!(processed, j)
            a, b = minmax(i, j)
            v = x[b] - x[a]
            nv = v / norm(v)
            d = dot(nv, (x[a] + x[b]) / 2)
            if i == a
                P, changed = clip(P, nv, d, j, tol)
            else
                P, changed = clip(P, -nv, -d, j, tol)
            end
            if changed
                isempty(P) && return P
                R = max_radius(P, xi)
            end
        end
        if finished || k >= n
            break
        end
        k = min(2 * k, n)
    end
    return P
end

# Local faces of a cell as (vertex loop, label)
_local_faces(P::ClipPolygon) = ((([k, k == length(P.verts) ? 1 : k + 1]), P.labels[k]) for k in eachindex(P.verts))
_local_faces(P::ClipPolyhedron) = zip(P.faces, P.labels)

function _assemble_cells(cells, merge_tol, ::Val{D}) where {D}
    nc = length(cells)
    pm = PointMerger{D}(merge_tol)
    faces = Vector{Vector{Int}}()
    neighbors = Tuple{Int, Int}[]
    boundary_tags = Int[]
    cell_face_lists = [Int[] for _ in 1:nc]
    creator = Dict{Tuple{Int, Int}, Int}()
    nmin = D == 2 ? 2 : 3
    for c in 1:nc
        P = cells[c]
        if isempty(P)
            @warn "Cell $c is empty"
            continue
        end
        g = [merge_point!(pm, SVector{D, Float64}(v)) for v in P.verts]
        for (loop, lab) in _local_faces(P)
            nodes = Int[]
            for k in loop
                id = g[k]
                if isempty(nodes) || nodes[end] != id
                    push!(nodes, id)
                end
            end
            while length(nodes) > 1 && nodes[end] == nodes[1]
                pop!(nodes)
            end
            length(unique(nodes)) < nmin && continue
            if lab < 0
                push!(faces, nodes)
                push!(neighbors, (c, 0))
                push!(boundary_tags, -lab)
                push!(cell_face_lists[c], length(faces))
            else
                j = lab
                key = minmax(c, j)
                cr = get(creator, key, 0)
                if cr != 0 && cr != c
                    # Face already created by the neighbor
                    continue
                end
                creator[key] = c
                if c < j
                    push!(faces, nodes)
                    push!(neighbors, (c, j))
                else
                    push!(faces, reverse(nodes))
                    push!(neighbors, (j, c))
                end
                push!(boundary_tags, 0)
                push!(cell_face_lists[c], length(faces))
                push!(cell_face_lists[j], length(faces))
            end
        end
    end
    # Compact nodes
    nn = length(pm.points)
    used = falses(nn)
    for f in faces, n in f
        used[n] = true
    end
    newid = cumsum(used)
    nodes = pm.points[used]
    faces_nodes = Int[]
    faces_nodespos = [1]
    for f in faces
        for n in f
            push!(faces_nodes, newid[n])
        end
        push!(faces_nodespos, length(faces_nodes) + 1)
    end
    cells_faces = Int[]
    cells_facepos = [1]
    for c in 1:nc
        append!(cells_faces, cell_face_lists[c])
        push!(cells_facepos, length(cells_faces) + 1)
    end
    return (nodes = nodes, faces_nodes = faces_nodes, faces_nodespos = faces_nodespos,
        neighbors = neighbors, cells_faces = cells_faces, cells_facepos = cells_facepos,
        boundary_tags = boundary_tags)
end
