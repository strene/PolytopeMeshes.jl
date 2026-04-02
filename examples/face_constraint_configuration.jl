# Face Constraint Configuration
#
# This example demonstrates how to configure the face constraint (fault)
# conformity method. It covers:
#   1. Interpolation vs exact representation of face constraints
#   2. Adapting the density in the constraint tessellation
#   3. Different fc_factor and circle_factor settings
#
# Based on: MRST book-ii/uprBookSection42.m (Section 4.2)

using PolytopeMeshes

## 1. Interpolation of face constraints
# Face constraints can be interpolated or represented exactly.
# Interpolation is useful when the path is a densely-sampled curve.
# Exact representation preserves "kinks" in the path.
println("=== 1. Interpolation vs Exact Representation ===")

# A curved constraint (100 points) and two L-shaped constraints
x_curve = collect(range(0.2, 0.8, length=100))
y_curve = 0.4 .+ 0.3 .* sin.(π .* x_curve)

lines = [Float64[x_curve y_curve],
         Float64[0.0 0.65; 0.0 0.4; 0.25 0.4],
         Float64[0.75 0.4; 1.0 0.4; 1.0 0.65]]

# Using surface_sites_2d directly to show the tessellation
F_interp = surface_sites_2d(lines, 0.1;
    interpolate_fc=[true, true, false])

println("  Curved constraint: $(length(x_curve)) input points")
println("  Interpolated → Circle centers: $(size(F_interp.c_CC, 1))")
println("  Face sites: $(size(F_interp.f_pts, 1))")
println("  Tip sites: $(size(F_interp.t_pts, 1))")
println()

## 2. Adaptive density along constraint
# A distance function can control the tessellation density.
println("=== 2. Adaptive Density Along Constraint ===")

F_adaptive = surface_sites_2d([lines[1]], 0.1;
    interpolate_fc=true,
    dist_fun=x -> 0.05 .+ 0.125 .* x[:, 1])

println("  Uniform density → Circle centers: $(size(F_interp.c_CC, 1))")
println("  Adaptive density → Circle centers: $(size(F_adaptive.c_CC, 1))")
println()

## 3. fc_factor comparison
# fc_factor controls the distance between face constraint sites relative
# to the background grid size.
println("=== 3. fc_factor Comparison ===")
f = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.5]]

G1, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, fc_factor=1.0)
G2, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, fc_factor=0.5)
G3, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, fc_factor=0.25)
println("  fc_factor=1.0  → Cells: $(G1.cells.num)")
println("  fc_factor=0.5  → Cells: $(G2.cells.num)")
println("  fc_factor=0.25 → Cells: $(G3.cells.num)")
println()

## 4. circle_factor comparison
# circle_factor controls the radius of the circles used to generate fault
# sites. Must be in (0.5, 1.0). Larger values create more elongated cells.
println("=== 4. circle_factor Comparison ===")

G1, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, circle_factor=0.55)
G2, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, circle_factor=0.7)
G3, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    face_constraints=f, circle_factor=0.9)
println("  circle_factor=0.55 → Cells: $(G1.cells.num)")
println("  circle_factor=0.7  → Cells: $(G2.cells.num)")
println("  circle_factor=0.9  → Cells: $(G3.cells.num)")
