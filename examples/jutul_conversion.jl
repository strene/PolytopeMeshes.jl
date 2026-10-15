# Convert a PolytopeMeshes grid to a Jutul UnstructuredMesh. The constraint
# tags become Jutul mesh entity tags.
using PolytopeMeshes
using Jutul

G = composite_pebi_grid_2d([0.05, 0.05], [1, 1];
    face_constraints = [[[0.2, 0.2], [0.8, 0.7]]],
    cell_constraints = [[[0.2, 0.8], [0.8, 0.2]]]
)
g = UnstructuredMesh(G)
geo = tpfv_geometry(g)
println("UnstructuredMesh with ", number_of_cells(g), " cells and ", number_of_faces(g), " interior faces")
println("Total area: ", sum(geo.volumes))
well_cells = Jutul.get_mesh_entity_tag(g, Cells(), :cell_constraints, :constraint_1)
fault_faces = Jutul.get_mesh_entity_tag(g, Faces(), :face_constraints, :constraint_1)
println("Well cells: ", well_cells)
println("Number of fault faces: ", length(fault_faces))
