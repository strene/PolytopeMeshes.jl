# Two intersecting faults and a deviated well in 3D (UPR example
# compositePebiGridWithTwoFaults). The faults are given as triangulated
# surfaces.
using PolytopeMeshes

nx = 10
x = range(0, 1, length = nx)
li = LinearIndices((nx, nx))
triangles = NTuple{3, Int}[]
for i in 1:(nx - 1), j in 1:(nx - 1)
    push!(triangles, (li[i, j], li[i + 1, j], li[i + 1, j + 1]))
    push!(triangles, (li[i, j], li[i + 1, j + 1], li[i, j + 1]))
end
fault1 = vec([[a, 0.5a + 0.1, b] .* 0.99 .+ 0.005 for a in x, b in x])
fault2 = vec([[a, -0.5a + 0.6, b] for a in x, b in x])
well = [[[0.9, 0.75, 0.1], [0.85, 0.85, 0.5], [0.9, 0.9, 0.9]]]

G = composite_voronoi_grid_3d((nx, nx, nx), (1, 1, 1);
    face_constraints = [(fault1, triangles), (fault2, triangles)],
    fc_rho = 1 / (nx - 1),
    cell_constraints = well,
    mlqt_max_level = 1
)
println(G)
println(check_mesh(G))
write_vtk("faults_and_well_3d", G)
write_vtk_faces("faults_and_well_3d_faults", G)

# A planar fault that ends inside the domain, triangulated automatically
S = surface_triangulation([[0.2, 0.2, 0.3], [0.8, 0.3, 0.4], [0.7, 0.8, 0.7], [0.2, 0.7, 0.6]], 0.08)
G = composite_voronoi_grid_3d((12, 12, 12), (1, 1, 1); face_constraints = [S], fc_rho = 0.08)
write_vtk("planar_fault_3d", G)
write_vtk_faces("planar_fault_3d_fault", G)
