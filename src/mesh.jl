"""
    CellData

Cell topology and geometry data.

# Fields
- `num`: Number of cells
- `faces`: Cell-to-face mapping (indices into face list)
- `facePos`: Position pointers for cell-face mapping; cell `i` has faces
  `faces[facePos[i]:facePos[i+1]-1]`
- `volumes`: Cell volumes (computed by `compute_geometry!`)
- `centroids`: Cell centroids (computed by `compute_geometry!`)
- `tag`: Boolean tags for constrained cells (e.g., well cells)
"""
mutable struct CellData
    num::Int
    faces::Vector{Int}
    facePos::Vector{Int}
    volumes::Vector{Float64}
    centroids::Matrix{Float64}
    tag::Vector{Bool}
end

"""
    FaceData

Face topology and geometry data.

# Fields
- `num`: Number of faces
- `nodes`: Face-to-node mapping (indices into node list)
- `nodePos`: Position pointers for face-node mapping; face `i` has nodes
  `nodes[nodePos[i]:nodePos[i+1]-1]`
- `neighbors`: Face neighbor cells [num_faces × 2]; 0 indicates boundary
- `areas`: Face areas (computed by `compute_geometry!`)
- `normals`: Face normal vectors (computed by `compute_geometry!`)
- `centroids`: Face centroids (computed by `compute_geometry!`)
- `tag`: Boolean tags for constrained faces (e.g., fault faces)
"""
mutable struct FaceData
    num::Int
    nodes::Vector{Int}
    nodePos::Vector{Int}
    neighbors::Matrix{Int}
    areas::Vector{Float64}
    normals::Matrix{Float64}
    centroids::Matrix{Float64}
    tag::Vector{Bool}
end

"""
    NodeData

Node coordinate data.

# Fields
- `num`: Number of nodes
- `coords`: Node coordinates [num_nodes × dim]
"""
mutable struct NodeData
    num::Int
    coords::Matrix{Float64}
end

"""
    UnstructuredMesh

An unstructured mesh data structure compatible with the MRST mesh format.

# Fields
- `cells`: Cell (volume element) information
- `faces`: Face (surface element) information
- `nodes`: Node (vertex) information
- `meshdim`: Mesh dimension (2 or 3)
- `type`: Mesh type identifier
"""
mutable struct UnstructuredMesh
    cells::CellData
    faces::FaceData
    nodes::NodeData
    meshdim::Int
    type::Vector{String}
end

"""
    UnstructuredMesh(nodes_coords, cell_faces, cell_facePos, face_nodes,
                     face_nodePos, face_neighbors; meshdim=2)

Construct an `UnstructuredMesh` from topology arrays.
"""
function UnstructuredMesh(
    node_coords::Matrix{Float64},
    cell_faces::Vector{Int},
    cell_facePos::Vector{Int},
    face_nodes::Vector{Int},
    face_nodePos::Vector{Int},
    face_neighbors::Matrix{Int};
    meshdim::Int = 2
)
    num_cells = length(cell_facePos) - 1
    num_faces = length(face_nodePos) - 1
    num_nodes = size(node_coords, 1)
    dim = size(node_coords, 2)

    cells = CellData(
        num_cells,
        cell_faces,
        cell_facePos,
        Float64[],
        zeros(0, dim),
        fill(false, num_cells)
    )

    faces = FaceData(
        num_faces,
        face_nodes,
        face_nodePos,
        face_neighbors,
        Float64[],
        zeros(0, dim),
        zeros(0, dim),
        fill(false, num_faces)
    )

    nodes = NodeData(num_nodes, node_coords)

    return UnstructuredMesh(cells, faces, nodes, meshdim, ["PolytopeMeshes"])
end

"""
    compute_geometry!(G::UnstructuredMesh)

Compute face areas, normals, centroids, and cell volumes and centroids.
Currently supports 2D meshes.
"""
function compute_geometry!(G::UnstructuredMesh)
    if G.meshdim == 2
        _compute_geometry_2d!(G)
    elseif G.meshdim == 3
        _compute_geometry_3d!(G)
    else
        error("Geometry computation not implemented for dimension $(G.meshdim)")
    end
    return G
end

function _compute_geometry_2d!(G::UnstructuredMesh)
    coords = G.nodes.coords
    nf = G.faces.num
    nc = G.cells.num

    # Face geometry
    face_centroids = zeros(nf, 2)
    face_normals = zeros(nf, 2)
    face_areas = zeros(nf)

    for i in 1:nf
        n1 = G.faces.nodes[G.faces.nodePos[i]]
        n2 = G.faces.nodes[G.faces.nodePos[i] + 1]
        p1 = coords[n1, :]
        p2 = coords[n2, :]

        face_centroids[i, :] = (p1 + p2) / 2
        d = p2 - p1
        face_normals[i, :] = [d[2], -d[1]]
        face_areas[i] = norm(d)
    end

    G.faces.areas = face_areas
    G.faces.normals = face_normals
    G.faces.centroids = face_centroids

    # Cell geometry (using the shoelace formula)
    cell_volumes = zeros(nc)
    cell_centroids = zeros(nc, 2)

    for i in 1:nc
        fi = G.cells.facePos[i]
        li = G.cells.facePos[i + 1] - 1

        # Collect ordered nodes for this cell
        cell_face_ids = G.cells.faces[fi:li]
        npts = length(cell_face_ids)

        # Get ordered nodes
        node_ids = Vector{Int}(undef, npts)
        for (j, fid) in enumerate(cell_face_ids)
            n1 = G.faces.nodes[G.faces.nodePos[fid]]
            n2 = G.faces.nodes[G.faces.nodePos[fid] + 1]
            # Choose node based on neighbor orientation
            if G.faces.neighbors[fid, 1] == i
                node_ids[j] = n1
            else
                node_ids[j] = n2
            end
        end

        # Shoelace formula
        area = 0.0
        cx = 0.0
        cy = 0.0
        for j in 1:npts
            k = mod(j, npts) + 1
            x1, y1 = coords[node_ids[j], 1], coords[node_ids[j], 2]
            x2, y2 = coords[node_ids[k], 1], coords[node_ids[k], 2]
            cross = x1 * y2 - x2 * y1
            area += cross
            cx += (x1 + x2) * cross
            cy += (y1 + y2) * cross
        end
        area /= 2.0
        cell_volumes[i] = abs(area)
        if abs(area) > eps()
            cell_centroids[i, :] = [cx, cy] / (6.0 * area)
        else
            # Degenerate cell - use centroid of nodes
            cell_centroids[i, :] = mean(coords[node_ids, :], dims=1)
        end
    end

    G.cells.volumes = cell_volumes
    G.cells.centroids = cell_centroids

    return G
end

function _compute_geometry_3d!(G::UnstructuredMesh)
    coords = G.nodes.coords
    nf = G.faces.num
    nc = G.cells.num

    # --- Face geometry ---
    face_centroids = zeros(nf, 3)
    face_normals = zeros(nf, 3)
    face_areas = zeros(nf)

    for i in 1:nf
        np_start = G.faces.nodePos[i]
        np_end = G.faces.nodePos[i + 1] - 1
        n_nodes = np_end - np_start + 1
        node_ids = G.faces.nodes[np_start:np_end]

        # Face centroid: area-weighted centroid of triangle fan
        # Face normal: sum of triangle cross products
        p1 = coords[node_ids[1], :]
        nx, ny, nz = 0.0, 0.0, 0.0
        cx, cy, cz = 0.0, 0.0, 0.0
        total_area = 0.0

        for j in 2:(n_nodes - 1)
            p2 = coords[node_ids[j], :]
            p3 = coords[node_ids[j + 1], :]
            # Cross product of triangle edges
            e1 = p2 - p1
            e2 = p3 - p1
            cross_x = e1[2] * e2[3] - e1[3] * e2[2]
            cross_y = e1[3] * e2[1] - e1[1] * e2[3]
            cross_z = e1[1] * e2[2] - e1[2] * e2[1]
            tri_area = 0.5 * sqrt(cross_x^2 + cross_y^2 + cross_z^2)
            # Triangle centroid
            tcx = (p1[1] + p2[1] + p3[1]) / 3.0
            tcy = (p1[2] + p2[2] + p3[2]) / 3.0
            tcz = (p1[3] + p2[3] + p3[3]) / 3.0
            cx += tri_area * tcx
            cy += tri_area * tcy
            cz += tri_area * tcz
            nx += cross_x
            ny += cross_y
            nz += cross_z
            total_area += tri_area
        end

        face_areas[i] = total_area
        if total_area > eps()
            face_centroids[i, :] = [cx, cy, cz] / total_area
        else
            # Degenerate face: use average of nodes
            face_centroids[i, :] = vec(mean(coords[node_ids, :], dims=1))
        end
        face_normals[i, :] = [nx, ny, nz] * 0.5
    end

    G.faces.areas = face_areas
    G.faces.normals = face_normals
    G.faces.centroids = face_centroids

    # --- Cell geometry ---
    # Use divergence theorem: V = (1/3) * sum_f sign_f * dot(normal_f, centroid_f)
    # Cell centroid via tetrahedral decomposition
    cell_volumes = zeros(nc)
    cell_centroids = zeros(nc, 3)

    for i in 1:nc
        fi = G.cells.facePos[i]
        li = G.cells.facePos[i + 1] - 1
        cell_face_ids = G.cells.faces[fi:li]

        # Reference point: average of face centroids
        p0 = zeros(3)
        for fid in cell_face_ids
            p0 .+= face_centroids[fid, :]
        end
        p0 ./= length(cell_face_ids)

        vol = 0.0
        cx, cy, cz = 0.0, 0.0, 0.0

        for fid in cell_face_ids
            # Determine sign: +1 if cell i is neighbor 1 (normal outward), -1 otherwise
            sign_f = G.faces.neighbors[fid, 1] == i ? 1.0 : -1.0

            np_start = G.faces.nodePos[fid]
            np_end = G.faces.nodePos[fid + 1] - 1
            node_ids = G.faces.nodes[np_start:np_end]
            n_nodes = length(node_ids)

            # Decompose face into triangles, form tetrahedra with p0
            v1 = coords[node_ids[1], :]
            for j in 2:(n_nodes - 1)
                v2 = coords[node_ids[j], :]
                v3 = coords[node_ids[j + 1], :]

                # Tetrahedron (p0, v1, v2, v3) - signed volume
                d1 = v1 - p0
                d2 = v2 - p0
                d3 = v3 - p0
                tet_vol = sign_f * (d1[1] * (d2[2] * d3[3] - d2[3] * d3[2]) -
                                    d1[2] * (d2[1] * d3[3] - d2[3] * d3[1]) +
                                    d1[3] * (d2[1] * d3[2] - d2[2] * d3[1])) / 6.0

                tet_c = (p0 + v1 + v2 + v3) / 4.0
                vol += tet_vol
                cx += tet_vol * tet_c[1]
                cy += tet_vol * tet_c[2]
                cz += tet_vol * tet_c[3]
            end
        end

        cell_volumes[i] = abs(vol)
        if abs(vol) > eps()
            cell_centroids[i, :] = [cx, cy, cz] / vol
        else
            cell_centroids[i, :] = p0
        end
    end

    G.cells.volumes = cell_volumes
    G.cells.centroids = cell_centroids

    return G
end
