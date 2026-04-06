module PolytopeMeshesGmshExt

using PolytopeMeshes
import Gmsh: gmsh

"""
    GmshTriangulation <: TriangulationBackend

Delaunay triangulation backend using [Gmsh](https://gmsh.info/).

This backend is automatically available when the `Gmsh` package is loaded
and is set as the default triangulation backend.

When used as a triangulation backend, Gmsh handles both the Delaunay
triangulation (via [`triangulate_2d`](@ref)) and full background grid
generation (via [`generate_background_grid_2d`](@ref)), bypassing
DistMesh entirely.

# Example
```julia
using PolytopeMeshes
using Gmsh  # activates GmshTriangulation as default

# Explicitly use the Gmsh backend:
G = clipped_voronoi_2d(pts, bnd; triangulation=GmshTriangulation())
```
"""
struct GmshTriangulation <: PolytopeMeshes.TriangulationBackend end

"""
    _gmsh_map_nodes_to_points(node_tags, node_coords, points)

Build a mapping from Gmsh node tags to original point indices by coordinate
proximity. Returns a `Dict{UInt64,Int}` mapping node tags to point indices.
"""
function _gmsh_map_nodes_to_points(node_tags, node_coords, points::Matrix{Float64})
    n = size(points, 1)
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
    return tag_to_idx
end

"""
    _gmsh_extract_triangles(elem_types, elem_node_tags, tag_to_output)

Extract triangles from Gmsh elements, remapping node tags to output point
indices using `tag_to_output`. Only triangles where all 3 nodes have a
mapping are included.

Returns an m×3 matrix of vertex indices.
"""
function _gmsh_extract_triangles(elem_types, elem_node_tags, tag_to_output::Dict)
    tri_idx = findfirst(==(2), elem_types)
    if tri_idx === nothing
        return zeros(Int, 0, 3)
    end

    tri_node_tags = elem_node_tags[tri_idx]
    n_tri = length(tri_node_tags) ÷ 3

    result = Vector{Vector{Int}}()
    for i in 1:n_tri
        t1 = tri_node_tags[3 * (i - 1) + 1]
        t2 = tri_node_tags[3 * (i - 1) + 2]
        t3 = tri_node_tags[3 * (i - 1) + 3]

        idx1 = get(tag_to_output, t1, 0)
        idx2 = get(tag_to_output, t2, 0)
        idx3 = get(tag_to_output, t3, 0)

        if idx1 > 0 && idx2 > 0 && idx3 > 0
            push!(result, [idx1, idx2, idx3])
        end
    end

    if isempty(result)
        return zeros(Int, 0, 3)
    end

    return reduce(vcat, [r' for r in result])
end

# ---------------------------------------------------------------------------
# triangulate_2d: pure Delaunay triangulation of an existing point set
# ---------------------------------------------------------------------------

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

        # Map Gmsh node tags to original point indices and extract triangles
        tag_to_idx = _gmsh_map_nodes_to_points(node_tags, node_coords, points)
        return _gmsh_extract_triangles(elem_types, elem_node_tags, tag_to_idx)
    finally
        gmsh.finalize()
    end
end

# ---------------------------------------------------------------------------
# generate_background_grid_2d: full mesh generation bypassing DistMesh
# ---------------------------------------------------------------------------

function PolytopeMeshes.generate_background_grid_2d(
    ::GmshTriangulation, fd::Function, fh::Function,
    h0::Float64, bbox::Matrix{Float64}, p_fix::Matrix{Float64};
    max_iter::Int = 500, poly_bdr::Matrix{Float64} = zeros(0, 2))

    nfix = size(p_fix, 1)

    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.model.add("background_grid_2d")

        # ── Create domain boundary geometry ──
        if size(poly_bdr, 1) >= 3
            # Polygon boundary
            n_bdr = size(poly_bdr, 1)
            bdr_tags = Vector{Int}(undef, n_bdr)
            for i in 1:n_bdr
                bdr_tags[i] = gmsh.model.geo.addPoint(
                    poly_bdr[i, 1], poly_bdr[i, 2], 0.0, h0)
            end
            lines = Vector{Int}(undef, n_bdr)
            for i in 1:n_bdr
                j = mod(i, n_bdr) + 1
                lines[i] = gmsh.model.geo.addLine(bdr_tags[i], bdr_tags[j])
            end
            cl = gmsh.model.geo.addCurveLoop(lines)
            s = gmsh.model.geo.addPlaneSurface([cl])
        else
            # Rectangular domain from bounding box
            x_min, y_min = bbox[1, 1], bbox[1, 2]
            x_max, y_max = bbox[2, 1], bbox[2, 2]
            p1 = gmsh.model.geo.addPoint(x_min, y_min, 0.0, h0)
            p2 = gmsh.model.geo.addPoint(x_max, y_min, 0.0, h0)
            p3 = gmsh.model.geo.addPoint(x_max, y_max, 0.0, h0)
            p4 = gmsh.model.geo.addPoint(x_min, y_max, 0.0, h0)
            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p3)
            l3 = gmsh.model.geo.addLine(p3, p4)
            l4 = gmsh.model.geo.addLine(p4, p1)
            cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            s = gmsh.model.geo.addPlaneSurface([cl])
        end

        # ── Embed fixed points with local mesh sizes ──
        if nfix > 0
            fix_tags = Vector{Int}(undef, nfix)
            for i in 1:nfix
                # Evaluate element size function at this point
                pt = reshape(p_fix[i:i, :], 1, 2)
                lc_i = h0 * max(fh(pt)[1], 0.01)
                fix_tags[i] = gmsh.model.geo.addPoint(
                    p_fix[i, 1], p_fix[i, 2], 0.0, lc_i)
            end
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.embed(0, fix_tags, 2, s)
        else
            gmsh.model.geo.synchronize()
        end

        # ── Configure meshing ──
        gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        # ── Generate mesh ──
        gmsh.model.mesh.generate(2)

        # ── Extract mesh data ──
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(2)

        if isempty(elem_types)
            return copy(p_fix), zeros(Int, 0, 3)
        end

        n_nodes = length(node_tags)

        # Build coordinate matrix and tag-to-position map
        all_pts = zeros(n_nodes, 2)
        tag_to_pos = Dict{UInt64,Int}()
        for i in eachindex(node_tags)
            tag_to_pos[node_tags[i]] = i
            all_pts[i, 1] = node_coords[3 * (i - 1) + 1]
            all_pts[i, 2] = node_coords[3 * (i - 1) + 2]
        end

        # Match fixed points to mesh nodes (by coordinate proximity)
        fix_node_pos = Dict{Int,Int}()  # mesh position → fixed point index
        for j in 1:nfix
            for i in 1:n_nodes
                if abs(all_pts[i, 1] - p_fix[j, 1]) < 1e-12 &&
                   abs(all_pts[i, 2] - p_fix[j, 2]) < 1e-12
                    fix_node_pos[i] = j
                    break
                end
            end
        end

        # ── Build output: fixed points first, then generated interior points ──
        # tag_to_output maps Gmsh node tags to output point indices
        tag_to_output = Dict{UInt64,Int}()
        # Tolerance for inside/outside tests — matches distmesh_2d's geps
        geps = 0.001 * h0

        # Map fixed point mesh nodes to their output indices
        for (mesh_pos, fix_idx) in fix_node_pos
            tag_to_output[node_tags[mesh_pos]] = fix_idx
        end

        # Collect generated interior points (not fixed, inside domain)
        gen_pts_list = Vector{Vector{Float64}}()
        next_idx = nfix + 1
        for i in 1:n_nodes
            if haskey(fix_node_pos, i)
                continue
            end
            pt = reshape(all_pts[i:i, :], 1, 2)
            d_val = fd(pt)[1]
            if d_val < geps
                push!(gen_pts_list, [all_pts[i, 1], all_pts[i, 2]])
                tag_to_output[node_tags[i]] = next_idx
                next_idx += 1
            end
        end

        # Assemble output point matrix
        if isempty(gen_pts_list)
            p_out = copy(p_fix)
        else
            gen_pts = reduce(vcat, [v' for v in gen_pts_list])
            p_out = vcat(p_fix, gen_pts)
        end

        # ── Extract triangles with remapped indices ──
        tri = _gmsh_extract_triangles(elem_types, elem_node_tags, tag_to_output)

        return p_out, tri
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
