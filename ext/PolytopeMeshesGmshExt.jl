module PolytopeMeshesGmshExt

using PolytopeMeshes
import Gmsh: gmsh

"""
    GmshTriangulation <: TriangulationBackend

Delaunay triangulation backend using [Gmsh](https://gmsh.info/).

This backend is automatically available when the `Gmsh` package is loaded
and is set as the default triangulation backend.

# Example
```julia
using PolytopeMeshes
using Gmsh  # activates GmshTriangulation as default

# Explicitly use the Gmsh backend:
G = clipped_voronoi_2d(pts, bnd; triangulation=GmshTriangulation())
```
"""
struct GmshTriangulation <: PolytopeMeshes.TriangulationBackend end

function PolytopeMeshes.triangulate_2d(::GmshTriangulation, points::Matrix{Float64})
    n = size(points, 1)
    if n < 3
        return zeros(Int, 0, 3)
    end

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.model.add("delaunay_2d")

        # Compute bounding box with large padding so that the bounding geometry
        # does not interfere with the Delaunay triangulation of the input points.
        min_x = minimum(points[:, 1])
        max_x = maximum(points[:, 1])
        min_y = minimum(points[:, 2])
        max_y = maximum(points[:, 2])
        d_max = max(max_x - min_x, max_y - min_y) + 1e-10
        pad = 100 * d_max
        lc = 1e6

        # Create bounding rectangle as geometry
        p1 = gmsh.model.geo.addPoint(min_x - pad, min_y - pad, 0.0, lc)
        p2 = gmsh.model.geo.addPoint(max_x + pad, min_y - pad, 0.0, lc)
        p3 = gmsh.model.geo.addPoint(max_x + pad, max_y + pad, 0.0, lc)
        p4 = gmsh.model.geo.addPoint(min_x - pad, max_y + pad, 0.0, lc)
        l1 = gmsh.model.geo.addLine(p1, p2)
        l2 = gmsh.model.geo.addLine(p2, p3)
        l3 = gmsh.model.geo.addLine(p3, p4)
        l4 = gmsh.model.geo.addLine(p4, p1)
        cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        s = gmsh.model.geo.addPlaneSurface([cl])

        # Add the input points to the geometry
        point_tags = Vector{Int}(undef, n)
        for i in 1:n
            point_tags[i] = gmsh.model.geo.addPoint(points[i, 1], points[i, 2], 0.0, lc)
        end

        gmsh.model.geo.synchronize()

        # Embed points in the surface so they become mesh nodes
        gmsh.model.mesh.embed(0, point_tags, 2, s)

        # Configure meshing: Delaunay algorithm, very large mesh size to
        # prevent Gmsh from inserting additional Steiner points.
        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(2)

        # Retrieve nodes and elements
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(2)

        if isempty(elem_types)
            return zeros(Int, 0, 3)
        end

        # Build mapping from Gmsh node tags to original point indices.
        # Match by coordinate proximity since Gmsh may reorder nodes.
        tag_to_idx = Dict{UInt64,Int}()
        for i in eachindex(node_tags)
            x = node_coords[3 * (i - 1) + 1]
            y = node_coords[3 * (i - 1) + 2]
            for j in 1:n
                if abs(x - points[j, 1]) < 1e-12 && abs(y - points[j, 2]) < 1e-12
                    tag_to_idx[node_tags[i]] = j
                    break
                end
            end
        end

        # Extract triangles (element type 2 = 3-node triangle)
        tri_idx = findfirst(==(2), elem_types)
        if tri_idx === nothing
            return zeros(Int, 0, 3)
        end

        tri_node_tags = elem_node_tags[tri_idx]
        n_tri = length(tri_node_tags) ÷ 3

        # Keep only triangles where all 3 nodes are original input points
        result = Vector{Vector{Int}}()
        for i in 1:n_tri
            t1 = tri_node_tags[3 * (i - 1) + 1]
            t2 = tri_node_tags[3 * (i - 1) + 2]
            t3 = tri_node_tags[3 * (i - 1) + 3]

            idx1 = get(tag_to_idx, t1, 0)
            idx2 = get(tag_to_idx, t2, 0)
            idx3 = get(tag_to_idx, t3, 0)

            if idx1 > 0 && idx2 > 0 && idx3 > 0
                push!(result, [idx1, idx2, idx3])
            end
        end

        if isempty(result)
            return zeros(Int, 0, 3)
        end

        return reduce(vcat, [r' for r in result])
    finally
        gmsh.finalize()
    end
end

# Export the new backend type
export GmshTriangulation

# Set Gmsh as the default backend when this extension loads
function __init__()
    PolytopeMeshes.set_default_triangulation_backend!(GmshTriangulation())
end

end # module PolytopeMeshesGmshExt
