# Cell Centroid Adaptation
#
# This example demonstrates how to adapt grids to cell-centroid constraints
# (wells). It covers:
#   1. Adaptation with and without protection layers
#   2. Combined face and cell constraints
#   3. Controlling MLQT refinement around wells
#
# Based on: MRST book-ii/uprBookSection5.m (Section 5)

using PolytopeMeshes
using LinearAlgebra

## 1. Protection layer comparison
# Protection sites are placed normal to the constraint path to create
# more rectangular cells around the well.
println("=== 1. Protection Layer Comparison ===")

constraint = [Float64[0.0 0.4; 0.2 0.5; 0.6 0.5; 0.8 0.6]]
grid_spacing = 0.12
distance_fn = p -> 0.12 * ones(size(p, 1))

# Without protection layer
cc_pts_no_prot, c_gs_no_prot, _, _ = line_sites_2d(constraint, grid_spacing;
    prot_layer=false)

# Create background sites
bg_x = collect(range(0.0, 1.0, length=10))
bg_y = collect(range(0.0, 1.0, length=10))
bg_sites = reduce(vcat, [[xi yi] for yi in bg_y for xi in bg_x])
bg_no_prot, _ = remove_conflict_points(bg_sites, cc_pts_no_prot, c_gs_no_prot)
sites_no_prot = vcat(cc_pts_no_prot, bg_no_prot)
bnd = Float64[0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
G_no_prot = clipped_pebi_2d(sites_no_prot, bnd)
compute_geometry!(G_no_prot)

# With protection layer
cc_pts_prot, c_gs_prot, prot_pts, p_gs = line_sites_2d(constraint, grid_spacing;
    prot_layer=true, prot_d=Function[distance_fn])
bg_prot, _ = remove_conflict_points(bg_sites, cc_pts_prot, c_gs_prot)
if size(prot_pts, 1) > 0
    bg_prot, _ = remove_conflict_points(bg_prot, prot_pts, p_gs)
end
sites_prot = vcat(cc_pts_prot, prot_pts, bg_prot)
G_prot = clipped_pebi_2d(sites_prot, bnd)
compute_geometry!(G_prot)

println("  Without protection:")
println("    Well sites: $(size(cc_pts_no_prot, 1))")
println("    Grid cells: $(G_no_prot.cells.num)")
println()
println("  With protection:")
println("    Well sites: $(size(cc_pts_prot, 1))")
println("    Protection sites: $(size(prot_pts, 1))")
println("    Grid cells: $(G_prot.cells.num)")
println()

## 2. Combined face and cell constraints
# Grid with both faults (face constraints) and wells (cell constraints).
println("=== 2. Combined Face and Cell Constraints ===")

face_lines = [Float64[0.2 0.2; 0.7 0.05],
              Float64[0.2 0.05; 0.7 0.2],
              Float64[0.1 0.4; 0.6 0.6],
              Float64[0.1 0.7; 0.45 0.7; 0.55 0.3]]
well_lines = [Float64[0.1 0.6; 0.2 0.6; 0.3 0.5; 0.4 0.3]]

G, pts, F = composite_pebi_grid_2d([0.05, 0.05], [1.0, 1.0];
    face_constraints=face_lines,
    cell_constraints=well_lines,
    cc_factor=0.5)
compute_geometry!(G)

println("  Face constraints: $(length(face_lines))")
println("  Cell constraints: $(length(well_lines))")
println("  Grid cells: $(G.cells.num)")
println("  Tagged cells (wells): $(count(G.cells.tag))")
println("  Tagged faces (faults): $(count(G.faces.tag))")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()

## 3. MLQT refinement control
# Different mlqt_max_level values control the refinement around wells.
println("=== 3. MLQT Refinement Levels ===")

w = [Float64[0.2 0.3; 0.5 0.5; 0.8 0.5]]
gS = [0.1, 0.1]
pdims = [1.0, 1.0]

configs = [
    (cc_factor=0.5,  mlqt=1, label="cc_factor=0.5,  mlqt=1"),
    (cc_factor=0.25, mlqt=2, label="cc_factor=0.25, mlqt=2"),
]

for cfg in configs
    local G
    G, _, _ = composite_pebi_grid_2d(gS, pdims;
        cell_constraints=w,
        cc_factor=cfg.cc_factor,
        mlqt_max_level=cfg.mlqt)
    compute_geometry!(G)
    println("  $(cfg.label) → Cells: $(G.cells.num), " *
            "Min vol: $(round(minimum(G.cells.volumes), sigdigits=3)), " *
            "Max vol: $(round(maximum(G.cells.volumes), sigdigits=3))")
end
