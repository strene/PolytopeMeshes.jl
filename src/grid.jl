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
    UnstructuredGrid

An unstructured grid data structure compatible with the MRST grid format.

# Fields
- `cells`: Cell (volume element) information
- `faces`: Face (surface element) information
- `nodes`: Node (vertex) information
- `griddim`: Grid dimension (2 or 3)
- `type`: Grid type identifier
"""
mutable struct UnstructuredGrid
    cells::CellData
    faces::FaceData
    nodes::NodeData
    griddim::Int
    type::Vector{String}
end

"""
    UnstructuredGrid(nodes_coords, cell_faces, cell_facePos, face_nodes,
                     face_nodePos, face_neighbors; griddim=2)

Construct an `UnstructuredGrid` from topology arrays.
"""
function UnstructuredGrid(
    node_coords::Matrix{Float64},
    cell_faces::Vector{Int},
    cell_facePos::Vector{Int},
    face_nodes::Vector{Int},
    face_nodePos::Vector{Int},
    face_neighbors::Matrix{Int};
    griddim::Int = 2
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

    return UnstructuredGrid(cells, faces, nodes, griddim, ["PolytopeMeshes"])
end

"""
    compute_geometry!(G::UnstructuredGrid)

Compute face areas, normals, centroids, and cell volumes and centroids.
Currently supports 2D grids.
"""
function compute_geometry!(G::UnstructuredGrid)
    if G.griddim == 2
        _compute_geometry_2d!(G)
    else
        error("3D geometry computation not yet implemented")
    end
    return G
end

function _compute_geometry_2d!(G::UnstructuredGrid)
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
