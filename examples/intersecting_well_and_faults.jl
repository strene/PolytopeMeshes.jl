# Intersecting Well and Faults
#
# This example demonstrates a single curved well path intersected by
# several straight faults. We use composite_pebi_grid_2d to create a PEBI
# grid conforming to both the well (cell constraint) and faults (face
# constraints).
#
# Based on: MRST 2d/intersectingWellAndFaults.m
# Reference: Berge et al. (2019), Computational Geosciences

using PolytopeMeshes

## Set well and fault paths
# The well is a curved path, and four straight faults cross it.

# Well: a parabolic curve y = 0.2 + x²
x = collect(range(0.1, 0.8, length=10))
y = 0.2 .+ x .^ 2
well = [Float64[x y]]

# Four faults crossing the well
fault1 = Float64[0.3 0.1; 0.1 0.3]
fault2 = Float64[0.3 0.1; 0.4 0.7]
fault3 = Float64[0.6 0.2; 0.3 0.7]
fault4 = Float64[0.8 0.4; 0.5 0.8]
faults = [fault1, fault2, fault3, fault4]

println("=== Intersecting Well and Faults ===")
println("  Well path: $(size(well[1], 1)) vertices")
println("  Number of faults: $(length(faults))")
println()

## Construct the grid
gS  = 1.0 / 24   # Grid size
wGf = 0.5         # Relative size of well cells compared to gS
fGf = 0.5         # Relative size of fault cells compared to gS
nRs = 1           # Number of MLQT refinement steps towards the well

G, pts, F = composite_pebi_grid_2d([gS, gS], [1.0, 1.0];
    cell_constraints=well,
    cc_factor=wGf,
    face_constraints=faults,
    fc_factor=fGf,
    mlqt_max_level=nRs)

compute_geometry!(G)

println("=== Grid Statistics ===")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Nodes: $(G.nodes.num)")
println("  Tagged cells (wells): $(count(G.cells.tag))")
println("  Tagged faces (faults): $(count(G.faces.tag))")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println("  Fault sites: $(size(F.f_pts, 1))")
