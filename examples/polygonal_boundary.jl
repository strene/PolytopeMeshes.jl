# Polygonal Boundary
#
# This example demonstrates creating PEBI grids on non-rectangular domains
# using the poly_bdr option of composite_pebi_grid_2d.
#
# Based on: Parts of MRST 2d/showOptionValuesCompositePebiGrid.m and
#           MRST book-ii/uprBookSection3.m

using PolytopeMeshes

## 1. Quadrilateral domain
# A non-axis-aligned quadrilateral boundary.
println("=== 1. Quadrilateral Domain ===")

bdr = Float64[0.0 0.0; -0.5 1.0; 0.5 1.5; 1.0 1.0]

G, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    poly_bdr=bdr)
compute_geometry!(G)

println("  Boundary vertices: $(size(bdr, 1))")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()

## 2. Star-shaped domain
# A more complex boundary shape.
println("=== 2. Star-shaped Domain ===")

bdr_star = Float64[
    0.0 0.0;
    0.5 0.2;
    1.0 0.0;
    0.8 0.5;
    1.0 1.0;
    0.5 0.8;
    0.0 1.0;
    0.2 0.5
]

G_star, _, _ = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    poly_bdr=bdr_star)
compute_geometry!(G_star)

println("  Boundary vertices: $(size(bdr_star, 1))")
println("  Cells: $(G_star.cells.num)")
println("  Faces: $(G_star.faces.num)")
println("  Total volume: $(round(sum(G_star.cells.volumes), digits=4))")
println()

## 3. Polygonal domain with faults
# Combine a non-rectangular domain with face constraints.
println("=== 3. Polygonal Domain with Faults ===")

bdr_poly = Float64[0.0 0.0; 1.0 0.0; 1.2 0.8; 0.8 1.2; -0.2 1.0]
faults = [Float64[0.2 0.3; 0.7 0.6],
          Float64[0.5 0.2; 0.3 0.8]]

G_poly, _, F = composite_pebi_grid_2d([0.08, 0.08], [1.0, 1.0];
    poly_bdr=bdr_poly,
    face_constraints=faults)
compute_geometry!(G_poly)

println("  Boundary vertices: $(size(bdr_poly, 1))")
println("  Face constraints: $(length(faults))")
println("  Cells: $(G_poly.cells.num)")
println("  Faces: $(G_poly.faces.num)")
println("  Tagged faces: $(count(G_poly.faces.tag))")
println("  Total volume: $(round(sum(G_poly.cells.volumes), digits=4))")
println()

## 4. Polygonal domain with wells
# Combine a non-rectangular domain with cell constraints.
println("=== 4. Polygonal Domain with Wells ===")

wells = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.7]]

G_well, _, _ = composite_pebi_grid_2d([0.08, 0.08], [1.0, 1.0];
    poly_bdr=bdr_poly,
    cell_constraints=wells,
    cc_factor=0.5)
compute_geometry!(G_well)

println("  Boundary vertices: $(size(bdr_poly, 1))")
println("  Cell constraints: $(length(wells))")
println("  Cells: $(G_well.cells.num)")
println("  Tagged cells (wells): $(count(G_well.cells.tag))")
println("  Total volume: $(round(sum(G_well.cells.volumes), digits=4))")
