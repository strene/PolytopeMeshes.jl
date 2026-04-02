# Clipped PEBI Grid from Arbitrary Sites
#
# This example demonstrates how to create a PEBI (Voronoi) grid directly
# from a set of generating sites, clipped to a polygonal boundary.
# This is the low-level interface; for most uses, composite_pebi_grid_2d
# is more convenient.
#
# Based on: MRST book-ii/uprBookSection2.m (Section 2.2, 2.3)

using PolytopeMeshes

## 1. Simple PEBI grid from a few sites
# Create a Voronoi tessellation from a small set of sites inside the
# unit square.
println("=== 1. Simple PEBI Grid ===")

sites = Float64[0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 1/3 1/3]
bnd = Float64[0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]

G = clipped_pebi_2d(sites, bnd)
compute_geometry!(G)

println("  Sites: $(size(sites, 1))")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Nodes: $(G.nodes.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()

## 2. PEBI grid from regular + perturbed sites
# Demonstrates the dual relationship between Delaunay triangulation
# and PEBI (Voronoi) grid.
println("=== 2. Regular + Perturbed Sites ===")

n = 5
sites_list = Float64[]
for j in 1:n
    for i in 1:n
        push!(sites_list, (i - 1) / (n - 1))
        push!(sites_list, (j - 1) / (n - 1))
    end
end
sites = collect(reshape(sites_list, 2, :)')

# Perturb interior points slightly
for i in 1:size(sites, 1)
    x, y = sites[i, :]
    if x > 0.01 && x < 0.99 && y > 0.01 && y < 0.99
        sites[i, :] .+= 0.05 .* randn(2)
    end
end

G = clipped_pebi_2d(sites, bnd)
compute_geometry!(G)

println("  Sites: $(size(sites, 1))")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()

## 3. Random sites
# Create a PEBI grid from random points.
println("=== 3. Random Sites ===")

using Random
Random.seed!(42)
n_random = 50
random_sites = rand(n_random, 2)
bnd_large = Float64[-0.1 -0.1; 1.1 -0.1; 1.1 1.1; -0.1 1.1]

G = clipped_pebi_2d(random_sites, bnd_large)
compute_geometry!(G)

println("  Random sites: $n_random")
println("  Cells: $(G.cells.num)")
println("  Faces: $(G.faces.num)")
println("  Total volume: $(round(sum(G.cells.volumes), digits=4))")
println()

# Verify all cells have positive volume
@assert all(G.cells.volumes .> 0) "All cells should have positive volume"
println("  ✓ All cells have positive volume")
println("  Min volume: $(round(minimum(G.cells.volumes), sigdigits=4))")
println("  Max volume: $(round(maximum(G.cells.volumes), sigdigits=4))")
