# Plotting of meshes with Makie. GLMakie gives interactive windows; replace it
# with CairoMakie to produce static figures (e.g. with `save("mesh.png", fig)`).
using PolytopeMeshes
using GLMakie

# 2D: a curved well crossing four faults
x = range(0.1, 0.8, length = 10)
well = [[[xi, 0.2 + xi^2] for xi in x]]
faults = [
    [[0.3, 0.1], [0.1, 0.3]],
    [[0.3, 0.1], [0.4, 0.7]],
    [[0.6, 0.2], [0.3, 0.7]],
    [[0.8, 0.4], [0.5, 0.8]],
]
G2 = pebi_grid_2d(1 / 24, [1, 1]; cell_constraints = well, cc_factor = 0.25,
    cc_refinement = true, cc_eps = 1 / 12, face_constraints = faults)

# Mesh with the faults (lines) and the well cells (markers) highlighted
fig, ax, plt = plot_polymesh(G2)
display(fig)
# Cell areas
fig, ax, plt = plot_polymesh_data(G2, cell_volumes(G2); colormap = :turbo)
display(fig)

# 3D: hierarchical gridding of three intersecting fractures
f1 = [[1, 3, 2], [4, 3, 2], [4, 3, 4], [1, 3, 4]]
f2 = [[2, 2, 3.3], [5, 2, 3.3], [5, 4, 3.3], [2, 4, 3.3]]
f3 = elliptic_surface_3d([3, 2.5, 3.2], 2, 2, 0, pi / 2, pi / 2)
fractures = [f1, f2, f3]
ints = surface_intersections_3d(fractures)
g2 = surface_grid_3d(fractures, line_grid_3d(ints, 0.2), ints, 0.25, 0.25 / 6)
G3 = volume_grid_3d([6, 6, 6], fractures, g2, 0.4, 0.4 / 12)

# Half of the cells colored by depth, with the fracture faces
fig, ax, plt = plot_polymesh(G3; cells = [s[2] < 3 for s in G3.cell_sites], color = [s[3] for s in G3.cell_sites])
display(fig)

# Only the fracture faces
fig = Figure()
ax = Axis3(fig[1, 1]; aspect = :data)
plot_polymesh_faces!(ax, G3)
display(fig)

# Interactive viewer: use the sliders to cut into the mesh, the menu to choose
# the cell field and the toggles to show the fractures without the cells.
plot_polymesh_interactive(G3; data = Dict("depth" => [s[3] for s in G3.cell_sites]))
plot_polymesh_interactive(G2)
