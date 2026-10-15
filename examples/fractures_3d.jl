# Hierarchical gridding of intersecting planar fractures (UPR example
# exact_pebi_gridding). The 1D grids of the fracture intersections are
# created first, then 2D grids of the fractures that conform to the 1D grids,
# and finally a 3D grid that conforms exactly to the 2D grids.
using PolytopeMeshes

f1 = [[1, 3, 2], [4, 3, 2], [4, 3, 4], [1, 3, 4]]
f2 = [[2, 2, 3.3], [5, 2, 3.3], [5, 4, 3.3], [2, 4, 3.3]]
f3 = elliptic_surface_3d([3, 2.5, 3.2], 2, 2, 0, pi / 2, pi / 2)
fractures = [f1, f2, f3]

intersections = surface_intersections_3d(fractures)
grids_1 = line_grid_3d(intersections, 0.2)
# 2D cell size and distance from the intersection lines to the 2D sites
ds_2 = 0.25
grids_2 = surface_grid_3d(fractures, grids_1, intersections, ds_2, ds_2 / 6)
# 3D cell size and distance from the fractures to the 3D sites
ds_3 = 0.4
G = volume_grid_3d([6, 6, 6], fractures, grids_2, ds_3, ds_3 / 12)
println(G)
println(check_mesh(G))
for k in eachindex(fractures)
    println("fracture $k: ", count(==(k), G.face_tags), " faces")
end
write_vtk("fractures_3d", G)
write_vtk_faces("fractures_3d_fractures", G)
