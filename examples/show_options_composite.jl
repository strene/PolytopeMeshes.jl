# Show Option Values for composite_voronoi_mesh_2d
#
# This example demonstrates all the major options for composite_voronoi_mesh_2d,
# corresponding to the MRST wrapper function compositePebiGrid2D.
#
# Based on: MRST 2d/showOptionValuesCompositePebiGrid.m

using PolytopeMeshes

## 1. cell_constraints
# Cell constraints are traced by cell centroids in the mesh.
# Each constraint is a piecewise-linear path (at least 2 vertices).
println("=== 1. Cell Constraints ===")
w = [Float64[0.2 0.8; 0.5 0.6; 0.8 0.8],
     Float64[0.5 0.2; 0.5 0.3]]  # short well segment
G, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w)
println("  Cells: $(G.cells.num), Tagged: $(count(G.cells.tag))")
println()

## 2. cc_factor
# Controls relative distance between sites along cell constraints.
# cc_factor=0.5 means half the reservoir mesh spacing.
println("=== 2. cc_factor ===")
w = [Float64[0.2 0.3; 0.8 0.7]]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, cc_factor=1.0)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, cc_factor=0.5)
println("  cc_factor=1.0 → Cells: $(G1.cells.num)")
println("  cc_factor=0.5 → Cells: $(G2.cells.num)")
println()

## 3. mlqt_max_level
# Controls number of MLQT (multi-level quadtree) refinement steps around
# cell constraints. More levels = smoother transition from fine to coarse.
println("=== 3. mlqt_max_level ===")
w = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.5]]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, cc_factor=0.5, mlqt_max_level=1)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, cc_factor=0.25, mlqt_max_level=2)
println("  mlqt_max_level=1 → Cells: $(G1.cells.num)")
println("  mlqt_max_level=2 → Cells: $(G2.cells.num)")
println()

## 4. cc_rho
# A function that controls the relative distance between cell constraint
# sites across the domain. cc_rho(p) should return a vector of relative
# distances for each point p.
println("=== 4. cc_rho ===")
w = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.7]]
cc_rho_custom = x -> 1.0 .- 0.9 .* x[:, 1]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, cc_rho=cc_rho_custom)
println("  Uniform cc_rho   → Cells: $(G1.cells.num)")
println("  Variable cc_rho  → Cells: $(G2.cells.num)")
println()

## 5. prot_layer
# Adds a protection layer of sites placed normal to the constraint path.
# This helps create more rectangular cells around the constraint.
println("=== 5. prot_layer ===")
x_vals = collect(range(0.2, 0.8, length=20))
y_vals = 0.5 .+ 0.1 .* sin.(π .* x_vals)
w = [Float64[x_vals y_vals]]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, interpolate_cc=true, prot_layer=false)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=w, interpolate_cc=true, prot_layer=true)
println("  Protection off → Cells: $(G1.cells.num)")
println("  Protection on  → Cells: $(G2.cells.num)")
println()

## 6. face_constraints
# Face constraints are traced by faces (edges) of the mesh.
# Used for faults, fractures, etc.
println("=== 6. face_constraints ===")
f = [Float64[0.2 0.8; 0.5 0.65; 0.8 0.8],
     Float64[0.5 0.2; 0.1 0.5]]

G, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f)
println("  Cells: $(G.cells.num), Tagged faces: $(count(G.faces.tag))")
println()

## 7. fc_factor
# Controls relative distance between face constraint sites.
println("=== 7. fc_factor ===")
f = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.5]]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, fc_factor=1.0)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, fc_factor=0.5)
println("  fc_factor=1.0 → Cells: $(G1.cells.num)")
println("  fc_factor=0.5 → Cells: $(G2.cells.num)")
println()

## 8. circle_factor
# Controls the radius of the circles used to generate face constraint sites.
# Must be in (0.5, 1.0). Larger values place sites further from the fault.
println("=== 8. circle_factor ===")
f = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.5]]

G1, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, circle_factor=0.55)
G2, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, circle_factor=0.9)
println("  circle_factor=0.55 → Cells: $(G1.cells.num)")
println("  circle_factor=0.9  → Cells: $(G2.cells.num)")
println()

## 9. poly_bdr
# Create a non-rectangular reservoir domain by specifying polygon vertices.
println("=== 9. poly_bdr (polygonal boundary) ===")
bdr = Float64[0.0 0.0; -0.5 1.0; 0.5 1.5; 1.0 1.0]

G, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
    poly_bdr=bdr)
compute_geometry!(G)
println("  Polygonal boundary with $(size(bdr, 1)) vertices")
println("  Cells: $(G.cells.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
