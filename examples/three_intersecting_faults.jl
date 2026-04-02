# Three Intersecting Faults
#
# This example shows three faults that intersect each other, inspired by
# the example in:
#   X. Y. Ding and L. S. K. Fung. "An unstructured Gridding Method for
#   Simulating Faulted Reservoirs populated with Complex Wells."
#   SPE Reservoir Simulation Symposium, Houston, Texas, February 2015.
#   doi: 10.2118/173243-MS
#
# Based on: MRST 2d/threeIntersectingFaults.m

using PolytopeMeshes

## Define the three fault paths
# Each fault is a piecewise-linear curve
fault1 = Float64[0.1 0.42; 0.4 0.55; 0.7 0.65]
fault2 = Float64[0.8 0.13; 0.6 0.4; 0.55 0.6]
fault3 = Float64[0.42 1.08; 0.45 0.9; 0.5 0.8; 0.58 0.6]

faults = [fault1, fault2, fault3]

println("=== Three Intersecting Faults ===")
for (i, f) in enumerate(faults)
    println("  Fault $i: $(size(f, 1)) vertices")
end
println()

## Create mesh
# Note: Use a moderate mesh size for reasonable computation time.
# For finer grids, consider using an external Delaunay package.
gS = 0.1  # Mesh size

G, pts, F = composite_voronoi_mesh_2d([gS, gS], [1.0, 1.15];
    face_constraints=faults)

compute_geometry!(G)

println("=== Mesh Statistics ===")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Nodes: $(G.nodes.num)")
println("  Tagged faces (faults): $(count(G.faces.tag))")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println("  Domain: [0, 1] × [0, 1.15]")
println()

# Verify the mesh tiles the domain correctly
expected_area = 1.0 * 1.15
actual_area = sum(G.cells.volumes)
println("  Expected area: $expected_area")
println("  Actual area:   $(round(actual_area, digits=4))")
println("  Area error:    $(round(abs(actual_area - expected_area) / expected_area * 100, digits=2))%")
