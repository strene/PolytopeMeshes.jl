# Well Branching
#
# This example demonstrates a branching well with several well
# intersections, using composite_pebi_grid_2d with cell constraints.
# Multi-level quadtree (MLQT) refinement is used to gradually refine
# the grid near the well paths.
#
# Based on: MRST 2d/wellBranching.m

using PolytopeMeshes

## Define well paths
# Each branch of the well is a separate cell constraint.
# Branches share common endpoints where they connect.

well1 = Float64[0.5 0.2; 0.5 0.3; 0.47 0.4; 0.4 0.5; 0.33 0.6; 0.26 0.7]
well2 = Float64[0.5 0.3; 0.53 0.4; 0.58 0.5]
well3 = Float64[0.5 0.45; 0.5 0.55; 0.45 0.65; 0.4 0.75; 0.38 0.85]
well4 = Float64[0.5 0.55; 0.55 0.65; 0.6 0.75; 0.62 0.85]

wells = [well1, well2, well3, well4]

println("=== Well Branching ===")
for (i, w) in enumerate(wells)
    println("  Branch $i: $(size(w, 1)) vertices")
end
println()

## Set gridding parameters
gS  = 1.0 / 24          # Grid size
wGf = 0.25               # Relative size of well cells
nRs = 2                  # Number of MLQT refinement levels

## Create composite PEBI grid
G, pts, F = composite_pebi_grid_2d([gS, gS], [1.0, 1.0];
    cell_constraints=wells,
    cc_factor=wGf,
    mlqt_max_level=nRs)

compute_geometry!(G)

println("=== Grid Statistics ===")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Nodes: $(G.nodes.num)")
println("  Tagged cells (wells): $(count(G.cells.tag))")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()
println("  Min cell volume: $(round(minimum(G.cells.volumes), sigdigits=4))")
println("  Max cell volume: $(round(maximum(G.cells.volumes), sigdigits=4))")
println("  Volume ratio:    $(round(maximum(G.cells.volumes) / minimum(G.cells.volumes), sigdigits=4))")
