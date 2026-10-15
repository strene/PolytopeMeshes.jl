# VTK output.

"""
    write_vtk(filename, mesh; cell_data = Dict())

Write the mesh to an unstructured VTK file (`.vtu`, appended to `filename`
if missing). 2D cells are written as polygons in the plane `z = 0`, 3D cells
as polyhedra. The cell tags, cell volumes and any extra `cell_data` (a
dictionary or named tuple of per-cell vectors) are written as cell data.
Returns the list of written files.
"""
function write_vtk(filename::AbstractString, m::PolytopalMesh{D}; cell_data = Dict{String, Any}()) where {D}
    points = zeros(3, number_of_nodes(m))
    for (i, p) in enumerate(m.node_points)
        points[1:D, i] = p
    end
    nc = number_of_cells(m)
    if D == 2
        cells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_POLYGON, cell_nodes(m, c)) for c in 1:nc]
    else
        cells = map(1:nc) do c
            faces = map(cell_faces(m, c)) do f
                nodes = collect(face_nodes(m, f))
                m.face_neighbors[f][1] == c ? Tuple(nodes) : Tuple(reverse(nodes))
            end
            WriteVTK.VTKPolyhedron(cell_nodes(m, c), faces...)
        end
    end
    geo = mesh_geometry(m)
    files = WriteVTK.vtk_grid(filename, points, cells) do vtk
        vtk["cell_tags", WriteVTK.VTKCellData()] = m.cell_tags
        vtk["volume", WriteVTK.VTKCellData()] = geo.cell_volumes
        vtk["cell_index", WriteVTK.VTKCellData()] = collect(1:nc)
        for (k, v) in pairs(cell_data)
            vtk[string(k), WriteVTK.VTKCellData()] = v
        end
    end
    return files
end

"""
    write_vtk_faces(filename, mesh; faces = findall(>(0), mesh.face_tags))

Write a subset of the faces (by default those on face constraints) to a VTK
file as polygons (3D) or lines (2D), with the face tags and boundary tags as
cell data. Useful to inspect the conformity to faults and fractures.
"""
function write_vtk_faces(filename::AbstractString, m::PolytopalMesh{D}; faces = findall(>(0), m.face_tags)) where {D}
    points = zeros(3, number_of_nodes(m))
    for (i, p) in enumerate(m.node_points)
        points[1:D, i] = p
    end
    ctype = D == 2 ? WriteVTK.VTKCellTypes.VTK_LINE : WriteVTK.VTKCellTypes.VTK_POLYGON
    cells = [WriteVTK.MeshCell(ctype, collect(face_nodes(m, f))) for f in faces]
    files = WriteVTK.vtk_grid(filename, points, cells) do vtk
        vtk["face_tags", WriteVTK.VTKCellData()] = m.face_tags[faces]
        vtk["boundary_tags", WriteVTK.VTKCellData()] = m.boundary_tags[faces]
        vtk["face_index", WriteVTK.VTKCellData()] = collect(faces)
    end
    return files
end
