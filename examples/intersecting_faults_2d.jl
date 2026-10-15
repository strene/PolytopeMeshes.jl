# Three intersecting faults inspired by Ding and Fung (SPE 173243) (UPR
# example threeIntersectingFaults), in a polygonal domain.
using PolytopeMeshes

faults = [
    [[0.1, 0.42], [0.4, 0.55], [0.7, 0.65]],
    [[0.8, 0.13], [0.6, 0.4], [0.55, 0.6]],
    [[0.42, 1.08], [0.45, 0.9], [0.5, 0.8], [0.58, 0.6]],
]
G = composite_pebi_grid_2d([1 / 30, 1 / 30], [1, 1.15]; face_constraints = faults)
println(G)
println(check_mesh(G))
write_vtk("intersecting_faults", G)
write_vtk_faces("intersecting_faults_faces", G)

# The same faults in a polygonal domain with a fully unstructured grid
polygon = [[0, 0], [1, 0], [1.1, 0.7], [0.6, 1.15], [0, 1]]
G = pebi_grid_2d(1 / 30, nothing; polygon = polygon, face_constraints = faults, fc_factor = 0.5)
println(G)
write_vtk("intersecting_faults_polygon", G)
