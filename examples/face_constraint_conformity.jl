# Face Constraint Conformity
#
# This example demonstrates how PEBI grids can conform to face constraints
# (faults/fractures). It shows both:
#   1. Manual construction using low-level functions
#   2. Automatic construction using composite_pebi_grid_2d
#
# Based on: MRST book-ii/uprBookSection41.m (Section 4.1)

using PolytopeMeshes
using LinearAlgebra

## Part 1: Making a constrained grid manually
# This demonstrates the basic principle of face-constraint conformity.
#
# We place circles at each vertex of the constraint path.
# The sites for the Voronoi grid are placed at the intersection points
# of consecutive circles, on both sides of the constraint. This ensures
# the resulting PEBI faces align with the constraint path.

println("=== Part 1: Manual Face Constraint Grid ===")

# Define a piecewise-linear face constraint
v = Float64[0.0 0.4; 0.2 0.5; 0.4 0.5; 0.6 0.6]

# Calculate distances between consecutive vertices
d = [norm(v[i+1, :] - v[i, :]) for i in 1:(size(v, 1) - 1)]
R = 0.6 * minimum(d)  # Circle radius

# Calculate circle intersections (face constraint sites)
# The sites are placed at the intersection of consecutive circles
dn = sqrt.(R^2 .- (d ./ 2) .^ 2)  # Normal offset

# Tangent and normal vectors along the constraint
t_vecs = zeros(length(d), 2)
n_vecs = zeros(length(d), 2)
for i in 1:length(d)
    t_vecs[i, :] = (v[i+1, :] - v[i, :]) / d[i]
    n_vecs[i, :] = [-t_vecs[i, 2], t_vecs[i, 1]]
end

# Sites at circle intersections
centers     = v[1:end-1, :] .+ (d ./ 2) .* t_vecs
left_sites  = centers .+ dn .* n_vecs
right_sites = centers .- dn .* n_vecs

# Tip site at the end of the constraint
tip_site = v[end, :] .+ R / sqrt(2)

# Combine all constraint sites
fault_sites = vcat(left_sites, right_sites, tip_site')

# Add background sites on a regular grid, removing conflicts
bg_sites = reduce(vcat, [[xi yi] for yi in 0.0:0.2:1.0 for xi in 0.0:0.2:1.0])

# Remove background sites that conflict with constraint sites
ref_pts = vcat(v, tip_site')
bg_sites, _ = remove_conflict_points(bg_sites, ref_pts,
    fill(R, size(ref_pts, 1)))

# Combine all sites
all_sites = vcat(fault_sites, bg_sites)

# Create clipped PEBI grid
bnd = Float64[0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
G = clipped_pebi_2d(all_sites, bnd)
compute_geometry!(G)

println("  Constraint vertices: $(size(v, 1))")
println("  Circle radius: $(round(R, digits=4))")
println("  Fault sites: $(size(fault_sites, 1))")
println("  Background sites: $(size(bg_sites, 1))")
println("  Total sites: $(size(all_sites, 1))")
println("  Grid cells: $(G.cells.num)")
println("  Grid faces: $(G.faces.num)")
println()

## Part 2: Automatic grid using composite_pebi_grid_2d
# The same kind of grid can be constructed automatically.
println("=== Part 2: Automatic Face Constraint Grid ===")

lines = [Float64[0.2 0.2; 0.7 0.05],
         Float64[0.2 0.05; 0.7 0.2],
         Float64[0.1 0.4; 0.6 0.6],
         Float64[0.1 0.7; 0.45 0.7; 0.55 0.3]]

G2, pts, F = composite_pebi_grid_2d([0.05, 0.05], [1.0, 1.0];
    face_constraints=lines)
compute_geometry!(G2)

println("  Face constraints: $(length(lines))")
println("  Grid cells: $(G2.cells.num)")
println("  Grid faces: $(G2.faces.num)")
println("  Tagged faces: $(count(G2.faces.tag))")
println("  Total volume: $(round(sum(G2.cells.volumes), digits=4))")
