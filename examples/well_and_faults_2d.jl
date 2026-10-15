# A curved well crossing four faults (UPR example intersectingWellAndFaults).
# The grid is created both with a Cartesian background grid
# (composite_voronoi_grid_2d) and a fully unstructured background grid
# (voronoi_grid_2d), and written to VTK for visualization in ParaView.
using PolytopeMeshes

x = range(0.1, 0.8, length = 10)
well = [[[xi, 0.2 + xi^2] for xi in x]]
faults = [
    [[0.3, 0.1], [0.1, 0.3]],
    [[0.3, 0.1], [0.4, 0.7]],
    [[0.6, 0.2], [0.3, 0.7]],
    [[0.8, 0.4], [0.5, 0.8]],
]
gs = 1 / 24 # Background grid size

# Cartesian background grid refined once around the well. Well and fault
# cells are half the size of the background cells.
Gc = composite_voronoi_grid_2d([gs, gs], [1, 1];
    cell_constraints = well, cc_factor = 0.5,
    face_constraints = faults, fc_factor = 0.5,
    mlqt_max_level = 1
)

# Unstructured background grid, refined towards the well.
Gd = voronoi_grid_2d(gs, [1, 1];
    cell_constraints = well, cc_factor = 0.25,
    cc_refinement = true, cc_eps = 1 / 12,
    face_constraints = faults
)

for (name, G) in (("well_and_faults_composite", Gc), ("well_and_faults_voronoi", Gd))
    println(name, ": ", G)
    # Cells along the well are tagged with the well index, faces on the
    # faults with the fault index.
    println("  well cells: ", count(>(0), G.cell_tags), ", fault faces: ", count(>(0), G.face_tags))
    write_vtk(name, G)
    write_vtk_faces(name * "_faults", G)
end
