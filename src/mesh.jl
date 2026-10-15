"""
    PolytopalMesh{D}

Polytopal mesh in `D` dimensions (`D = 2` or `3`) stored in a face-based
format compatible with MRST and `Jutul.UnstructuredMesh`.

# Fields
- `node_points`: node coordinates.
- `faces_nodes`, `faces_nodespos`: nodes of face `f` are
  `faces_nodes[faces_nodespos[f]:faces_nodespos[f+1]-1]`. In 2D a face is an
  edge with two nodes. In 3D the nodes are ordered so that the right-hand
  normal points from the first to the second neighbor (out of the domain for
  boundary faces). In 2D the normal of the edge `a → b` is `(b - a)` rotated
  clockwise.
- `face_neighbors`: `(l, r)` cells of each face, with `r == 0` on the boundary.
- `cells_faces`, `cells_facepos`: faces of each cell.
- `cell_sites`: the Voronoi site (generating point) of each cell.
- `cell_tags`: index of the cell constraint (line/point) the cell was generated
  from, `0` otherwise.
- `face_tags`: index of the face constraint (surface) the face conforms to,
  `0` otherwise.
- `boundary_tags`: index of the domain facet a boundary face lies on, `0` for
  interior faces.
"""
struct PolytopalMesh{D}
    node_points::Vector{SVector{D, Float64}}
    faces_nodes::Vector{Int}
    faces_nodespos::Vector{Int}
    face_neighbors::Vector{Tuple{Int, Int}}
    cells_faces::Vector{Int}
    cells_facepos::Vector{Int}
    cell_sites::Vector{SVector{D, Float64}}
    cell_tags::Vector{Int}
    face_tags::Vector{Int}
    boundary_tags::Vector{Int}
end

number_of_cells(m::PolytopalMesh) = length(m.cells_facepos) - 1
number_of_faces(m::PolytopalMesh) = length(m.faces_nodespos) - 1
number_of_nodes(m::PolytopalMesh) = length(m.node_points)

face_nodes(m::PolytopalMesh, f) = view(m.faces_nodes, m.faces_nodespos[f]:(m.faces_nodespos[f + 1] - 1))
cell_faces(m::PolytopalMesh, c) = view(m.cells_faces, m.cells_facepos[c]:(m.cells_facepos[c + 1] - 1))
is_boundary_face(m::PolytopalMesh, f) = m.face_neighbors[f][2] == 0

function Base.show(io::IO, m::PolytopalMesh{D}) where {D}
    nb = count(x -> x[2] == 0, m.face_neighbors)
    print(io, "PolytopalMesh{$D} with $(number_of_cells(m)) cells, ")
    print(io, "$(number_of_faces(m)) faces ($nb boundary) and $(number_of_nodes(m)) nodes")
end

"""
    cell_nodes(m, c)

Nodes of cell `c`. In 2D the nodes are returned as a counter-clockwise loop,
in 3D as a sorted list of unique nodes.
"""
function cell_nodes(m::PolytopalMesh{2}, c)
    # Collect the edges oriented counter-clockwise around the cell and walk them.
    next = Dict{Int, Int}()
    for f in cell_faces(m, c)
        a, b = face_nodes(m, f)
        if m.face_neighbors[f][1] != c
            a, b = b, a
        end
        next[a] = b
    end
    isempty(next) && return Int[]
    start = first(keys(next))
    loop = [start]
    cur = next[start]
    while cur != start && length(loop) <= length(next)
        push!(loop, cur)
        cur = get(next, cur, start)
    end
    return loop
end

function cell_nodes(m::PolytopalMesh{3}, c)
    nodes = Int[]
    for f in cell_faces(m, c)
        append!(nodes, face_nodes(m, f))
    end
    return sort!(unique!(nodes))
end

# Face geometry: (centroid, measure, unit normal)
function _face_geometry(m::PolytopalMesh{2}, f)
    a, b = face_nodes(m, f)
    pa, pb = m.node_points[a], m.node_points[b]
    v = pb - pa
    L = norm(v)
    n = SVector(v[2], -v[1]) / L
    return (pa + pb) / 2, L, n
end

function _face_geometry(m::PolytopalMesh{3}, f)
    nodes = face_nodes(m, f)
    pts = m.node_points
    center = zero(Vec3)
    for n in nodes
        center += pts[n]
    end
    center /= length(nodes)
    area = 0.0
    centroid = zero(Vec3)
    normal = zero(Vec3)
    k = length(nodes)
    for i in 1:k
        a = pts[nodes[i]]
        b = pts[nodes[mod1(i + 1, k)]]
        cr = cross(a - center, b - center)
        A = norm(cr) / 2
        normal += cr
        area += A
        centroid += A * (a + b + center) / 3
    end
    if area > 0
        centroid /= area
    else
        centroid = center
    end
    return centroid, area, normal / max(norm(normal), eps())
end

"""
    mesh_geometry(m)

Compute geometric quantities of the mesh. Returns a named tuple with fields
`cell_volumes`, `cell_centroids`, `face_areas`, `face_normals` and
`face_centroids`. In 2D, "volumes" are areas and "areas" are edge lengths.
"""
function mesh_geometry(m::PolytopalMesh{D}) where {D}
    nf = number_of_faces(m)
    nc = number_of_cells(m)
    T = SVector{D, Float64}
    face_centroids = Vector{T}(undef, nf)
    face_areas = Vector{Float64}(undef, nf)
    face_normals = Vector{T}(undef, nf)
    for f in 1:nf
        c, a, n = _face_geometry(m, f)
        face_centroids[f] = c
        face_areas[f] = a
        face_normals[f] = n
    end
    cell_volumes = zeros(nc)
    cell_centroids = Vector{T}(undef, nc)
    for c in 1:nc
        # Approximate center used as apex of the sub-simplices
        x0 = zero(T)
        nfc = 0
        for f in cell_faces(m, c)
            x0 += face_centroids[f]
            nfc += 1
        end
        x0 /= max(nfc, 1)
        vol = 0.0
        cent = zero(T)
        for f in cell_faces(m, c)
            sgn = m.face_neighbors[f][1] == c ? 1.0 : -1.0
            nodes = face_nodes(m, f)
            if D == 2
                a = m.node_points[nodes[1]]
                b = m.node_points[nodes[2]]
                # Outward normal of a -> b is to the right, i.e. the edge is
                # traversed counter-clockwise around the cell.
                v = sgn * ((a[1] - x0[1]) * (b[2] - x0[2]) - (a[2] - x0[2]) * (b[1] - x0[1])) / 2
                vol += v
                cent += v * (x0 + a + b) / 3
            else
                xf = face_centroids[f]
                k = length(nodes)
                for i in 1:k
                    a = m.node_points[nodes[i]]
                    b = m.node_points[nodes[mod1(i + 1, k)]]
                    v = sgn * dot(cross(a - xf, b - xf), xf - x0) / 6
                    vol += v
                    cent += v * (x0 + xf + a + b) / 4
                end
            end
        end
        cell_volumes[c] = vol
        cell_centroids[c] = vol != 0 ? cent / vol : x0
    end
    return (
        cell_volumes = cell_volumes,
        cell_centroids = cell_centroids,
        face_areas = face_areas,
        face_normals = face_normals,
        face_centroids = face_centroids,
    )
end

cell_volumes(m::PolytopalMesh) = mesh_geometry(m).cell_volumes
cell_centroids(m::PolytopalMesh) = mesh_geometry(m).cell_centroids
face_areas(m::PolytopalMesh) = mesh_geometry(m).face_areas
face_normals(m::PolytopalMesh) = mesh_geometry(m).face_normals
face_centroids(m::PolytopalMesh) = mesh_geometry(m).face_centroids

"""
    check_mesh(m; rtol = 1e-8)

Check the consistency of a mesh: positive cell volumes, closed cells (the
oriented face normals of each cell sum to zero) and valid face neighbors.
Returns a named tuple with fields `ok`, `negative_volumes`, `open_cells` and
`max_closure_error`.
"""
function check_mesh(m::PolytopalMesh{D}; rtol = 1e-8) where {D}
    geo = mesh_geometry(m)
    nc = number_of_cells(m)
    open_cells = Int[]
    negative = findall(v -> v <= 0, geo.cell_volumes)
    max_err = 0.0
    for c in 1:nc
        s = zero(SVector{D, Float64})
        tot = 0.0
        for f in cell_faces(m, c)
            sgn = m.face_neighbors[f][1] == c ? 1.0 : -1.0
            s += sgn * geo.face_areas[f] * geo.face_normals[f]
            tot += geo.face_areas[f]
        end
        err = norm(s) / max(tot, eps())
        max_err = max(max_err, err)
        if err > rtol
            push!(open_cells, c)
        end
    end
    nb_ok = all(((l, r),) -> 1 <= l <= nc && 0 <= r <= nc && l != r, m.face_neighbors)
    ok = isempty(open_cells) && isempty(negative) && nb_ok
    return (ok = ok, negative_volumes = negative, open_cells = open_cells, max_closure_error = max_err)
end
