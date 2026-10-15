module PolytopeMeshesJutulExt

using PolytopeMeshes
using PolytopeMeshes: PolytopalMesh, number_of_cells, number_of_faces, face_nodes, cell_faces
import Jutul

Jutul.number_of_cells(m::PolytopalMesh) = number_of_cells(m)
Jutul.number_of_faces(m::PolytopalMesh) = number_of_faces(m)

"""
    Jutul.UnstructuredMesh(m::PolytopalMesh; tags = true, kwarg...)

Convert a [`PolytopalMesh`](@ref) to a `Jutul.UnstructuredMesh`. If `tags` is
true, the constraint tags are added as mesh entity tags:
- cells in group `:cell_constraints` with values `:constraint_<k>`,
- interior faces in group `:face_constraints` with values `:constraint_<k>`,
- boundary faces in group `:boundary` with values `:boundary_<k>`.
"""
function Jutul.UnstructuredMesh(m::PolytopalMesh{D}; tags = true, kwarg...) where {D}
    nc = number_of_cells(m)
    nf = number_of_faces(m)
    interior = findall(f -> m.face_neighbors[f][2] != 0, 1:nf)
    boundary = findall(f -> m.face_neighbors[f][2] == 0, 1:nf)
    int_index = zeros(Int, nf)
    bnd_index = zeros(Int, nf)
    int_index[interior] = 1:length(interior)
    bnd_index[boundary] = 1:length(boundary)

    function pack(faces)
        nodes = Int[]
        pos = [1]
        for f in faces
            append!(nodes, face_nodes(m, f))
            push!(pos, length(nodes) + 1)
        end
        return nodes, pos
    end
    faces_nodes, faces_nodespos = pack(interior)
    bnd_faces_nodes, bnd_faces_nodespos = pack(boundary)

    cells_faces = Int[]
    cells_facepos = [1]
    bnd_cells_faces = Int[]
    bnd_cells_facepos = [1]
    for c in 1:nc
        for f in cell_faces(m, c)
            if int_index[f] > 0
                push!(cells_faces, int_index[f])
            else
                push!(bnd_cells_faces, bnd_index[f])
            end
        end
        push!(cells_facepos, length(cells_faces) + 1)
        push!(bnd_cells_facepos, length(bnd_cells_faces) + 1)
    end
    neighbors = m.face_neighbors[interior]
    bnd_cells = [m.face_neighbors[f][1] for f in boundary]
    g = Jutul.UnstructuredMesh(
        cells_faces, cells_facepos,
        bnd_cells_faces, bnd_cells_facepos,
        faces_nodes, faces_nodespos,
        bnd_faces_nodes, bnd_faces_nodespos,
        copy(m.node_points),
        neighbors,
        bnd_cells;
        kwarg...
    )
    if tags
        for k in sort!(unique(m.cell_tags))
            k > 0 || continue
            Jutul.set_mesh_entity_tag!(g, Jutul.Cells(), :cell_constraints, Symbol("constraint_$k"), findall(==(k), m.cell_tags))
        end
        ftags = m.face_tags[interior]
        for k in sort!(unique(ftags))
            k > 0 || continue
            Jutul.set_mesh_entity_tag!(g, Jutul.Faces(), :face_constraints, Symbol("constraint_$k"), findall(==(k), ftags))
        end
        btags = m.boundary_tags[boundary]
        for k in sort!(unique(btags))
            k > 0 || continue
            Jutul.set_mesh_entity_tag!(g, Jutul.BoundaryFaces(), :boundary, Symbol("boundary_$k"), findall(==(k), btags))
        end
    end
    return g
end

end # module
