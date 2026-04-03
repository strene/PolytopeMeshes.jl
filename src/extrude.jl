"""
    extrude(G::UnstructuredMesh, length::Real, layers::Int; dim=:z)

Extrude a 2D `UnstructuredMesh` into 3D by creating `layers` uniform layers
up to the given `length` in the extrusion direction.

# Arguments
- `G`: A 2D `UnstructuredMesh` to extrude
- `length`: The total extent in the extrusion direction
- `layers`: Number of layers in the extrusion

# Keyword Arguments
- `dim=:z`: Extrusion direction (`:x`, `:y`, or `:z`). When `:x` or `:y`,
  the 2D coordinates are placed in the yz or xz plane respectively before
  extruding.

# Returns
A 3D `UnstructuredMesh`.
"""
function extrude(G::UnstructuredMesh, length::Real, layers::Int; dim::Symbol=:z)
    extrude(G, [Float64(length)], [layers]; dim=dim)
end

"""
    extrude(G::UnstructuredMesh, lengths::Vector{<:Real}, layers::Vector{<:Integer}; dim=:z)

Extrude a 2D `UnstructuredMesh` into 3D with variable layer spacing.

# Arguments
- `G`: A 2D `UnstructuredMesh` to extrude
- `lengths`: Coordinates of the mesh layers in the extrusion dimension. Each
  entry is the end coordinate of a segment.
- `layers`: Number of layers between each extrusion length. Must have the
  same length as `lengths`.

# Keyword Arguments
- `dim=:z`: Extrusion direction (`:x`, `:y`, or `:z`). When `:x` or `:y`,
  the 2D coordinates are placed in the yz or xz plane respectively before
  extruding.

# Returns
A 3D `UnstructuredMesh`.

# Example
```julia
# Extrude a 2D mesh to 3D with 2 layers from 0 to 1 and 3 layers from 1 to 3
G3D = extrude(G2D, [1.0, 3.0], [2, 3])
```
"""
function extrude(G::UnstructuredMesh, lengths::Vector{<:Real}, layers::Vector{<:Integer}; dim::Symbol=:z)
    @assert G.meshdim == 2 "Can only extrude 2D meshes"
    @assert length(lengths) == length(layers) "lengths and layers must have the same length"
    @assert all(l -> l >= 1, layers) "All layer counts must be at least 1"
    @assert dim in (:x, :y, :z) "dim must be :x, :y, or :z"

    # Compute z-coordinates for all levels
    z = Float64[0.0]
    prev_z = 0.0
    for i in eachindex(lengths)
        target = Float64(lengths[i])
        n = layers[i]
        for j in 1:n
            push!(z, prev_z + (target - prev_z) * j / n)
        end
        prev_z = target
    end

    nlevels = length(z)     # number of z-levels
    nlayers = nlevels - 1   # number of layers

    nn2d = G.nodes.num
    nf2d = G.faces.num
    nc2d = G.cells.num

    # --- 3D Nodes ---
    nn3d = nn2d * nlevels
    coords3d = zeros(nn3d, 3)
    for k in 1:nlevels
        for i in 1:nn2d
            idx = (k - 1) * nn2d + i
            coords3d[idx, 1] = G.nodes.coords[i, 1]
            coords3d[idx, 2] = G.nodes.coords[i, 2]
            coords3d[idx, 3] = z[k]
        end
    end

    # Apply coordinate transformation for dim
    if dim == :x
        # 2D (col1, col2) → 3D (extrude, col1, col2)
        coords3d = coords3d[:, [3, 1, 2]]
    elseif dim == :y
        # 2D (col1, col2) → 3D (col1, extrude, col2)
        coords3d = coords3d[:, [1, 3, 2]]
    end
    # dim == :z → no transformation needed

    # --- Precompute ordered cell nodes for each 2D cell ---
    cell_nodes_2d = Vector{Vector{Int}}(undef, nc2d)
    for c in 1:nc2d
        fi = G.cells.facePos[c]
        li = G.cells.facePos[c + 1] - 1
        cell_face_ids = G.cells.faces[fi:li]
        npts = length(cell_face_ids)
        node_ids = Vector{Int}(undef, npts)
        for (j, fid) in enumerate(cell_face_ids)
            n1 = G.faces.nodes[G.faces.nodePos[fid]]
            n2 = G.faces.nodes[G.faces.nodePos[fid] + 1]
            if G.faces.neighbors[fid, 1] == c
                node_ids[j] = n1
            else
                node_ids[j] = n2
            end
        end
        cell_nodes_2d[c] = node_ids
    end

    # --- 3D Faces ---
    # Type 1: Horizontal faces (polygons) - one per 2D cell per z-level
    nf_horiz = nc2d * nlevels
    # Type 2: Vertical faces (quads) - one per 2D face per layer
    nf_vert = nf2d * nlayers
    nf3d = nf_horiz + nf_vert

    # Build face-node connectivity
    face_nodes_vec = Int[]
    face_nodePos_vec = Int[1]

    # Horizontal faces: index (k-1)*nc2d + c for level k, cell c
    for k in 1:nlevels
        for c in 1:nc2d
            cnodes = cell_nodes_2d[c]
            for ni in cnodes
                push!(face_nodes_vec, (k - 1) * nn2d + ni)
            end
            push!(face_nodePos_vec, face_nodePos_vec[end] + length(cnodes))
        end
    end

    # Vertical faces: index nf_horiz + (k-1)*nf2d + f for layer k, face f
    for k in 1:nlayers
        for f in 1:nf2d
            n1 = G.faces.nodes[G.faces.nodePos[f]]
            n2 = G.faces.nodes[G.faces.nodePos[f] + 1]
            # Quad: [n1_k, n2_k, n2_(k+1), n1_(k+1)]
            push!(face_nodes_vec, (k - 1) * nn2d + n1)
            push!(face_nodes_vec, (k - 1) * nn2d + n2)
            push!(face_nodes_vec, k * nn2d + n2)
            push!(face_nodes_vec, k * nn2d + n1)
            push!(face_nodePos_vec, face_nodePos_vec[end] + 4)
        end
    end

    # --- Face Neighbors ---
    face_neighbors = zeros(Int, nf3d, 2)

    # Horizontal faces: neighbor 1 = cell below, neighbor 2 = cell above
    for k in 1:nlevels
        for c in 1:nc2d
            fidx = (k - 1) * nc2d + c
            if k == 1
                # Bottom boundary
                face_neighbors[fidx, 1] = 0
                face_neighbors[fidx, 2] = c  # cell (c, 1)
            elseif k == nlevels
                # Top boundary
                face_neighbors[fidx, 1] = (k - 2) * nc2d + c  # cell (c, nlayers)
                face_neighbors[fidx, 2] = 0
            else
                face_neighbors[fidx, 1] = (k - 2) * nc2d + c  # cell (c, k-1)
                face_neighbors[fidx, 2] = (k - 1) * nc2d + c  # cell (c, k)
            end
        end
    end

    # Vertical faces: inherit 2D neighbor structure
    for k in 1:nlayers
        for f in 1:nf2d
            fidx = nf_horiz + (k - 1) * nf2d + f
            c1_2d = G.faces.neighbors[f, 1]
            c2_2d = G.faces.neighbors[f, 2]
            face_neighbors[fidx, 1] = c1_2d > 0 ? (k - 1) * nc2d + c1_2d : 0
            face_neighbors[fidx, 2] = c2_2d > 0 ? (k - 1) * nc2d + c2_2d : 0
        end
    end

    # --- 3D Cells ---
    nc3d = nc2d * nlayers

    cell_faces_vec = Int[]
    cell_facePos_vec = Int[1]

    for k in 1:nlayers
        for c in 1:nc2d
            # Bottom horizontal face
            push!(cell_faces_vec, (k - 1) * nc2d + c)
            # Top horizontal face
            push!(cell_faces_vec, k * nc2d + c)
            # Side vertical faces
            fi = G.cells.facePos[c]
            li = G.cells.facePos[c + 1] - 1
            nside = li - fi + 1
            for fp in fi:li
                f2d = G.cells.faces[fp]
                push!(cell_faces_vec, nf_horiz + (k - 1) * nf2d + f2d)
            end
            push!(cell_facePos_vec, cell_facePos_vec[end] + 2 + nside)
        end
    end

    # --- Build face tags ---
    face_tags = fill(false, nf3d)
    # Propagate 2D face tags to vertical faces
    for k in 1:nlayers
        for f in 1:nf2d
            if G.faces.tag[f]
                face_tags[nf_horiz + (k - 1) * nf2d + f] = true
            end
        end
    end

    # --- Build cell tags ---
    cell_tags = fill(false, nc3d)
    # Propagate 2D cell tags to all layers
    for k in 1:nlayers
        for c in 1:nc2d
            if G.cells.tag[c]
                cell_tags[(k - 1) * nc2d + c] = true
            end
        end
    end

    # --- Construct the 3D mesh ---
    G3d = UnstructuredMesh(
        coords3d,
        cell_faces_vec,
        cell_facePos_vec,
        face_nodes_vec,
        face_nodePos_vec,
        face_neighbors;
        meshdim=3
    )

    G3d.faces.tag = face_tags
    G3d.cells.tag = cell_tags

    return G3d
end
