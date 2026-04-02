module PolytopeMeshesGLMakieExt

using PolytopeMeshes
using GLMakie

"""
    plot_mesh(G::UnstructuredMesh; kwargs...)

Plot a 2D `UnstructuredMesh` using GLMakie.

# Keyword Arguments
- `color_cells=true`: Color cells by index when `true`
- `colormap=:viridis`: Colormap for cell coloring
- `strokecolor=:black`: Color of cell edges
- `strokewidth=1.0`: Width of cell edges
- `show_nodes=false`: Show mesh nodes as scatter points
- `node_color=:red`: Color of node markers
- `node_size=4`: Size of node markers
- `figure_size=(800, 600)`: Size of the figure
- `title=""`: Title for the plot

Returns a `(fig, ax)` tuple.
"""
function PolytopeMeshes.plot_mesh(
    G::UnstructuredMesh;
    color_cells::Bool = true,
    colormap = :viridis,
    strokecolor = :black,
    strokewidth::Real = 1.0,
    show_nodes::Bool = false,
    node_color = :red,
    node_size::Real = 4,
    figure_size::Tuple{Int,Int} = (800, 600),
    title::String = "",
)
    coords = G.nodes.coords
    nc = G.cells.num

    fig = Figure(; size = figure_size)
    ax = Axis(fig[1, 1]; title = title, aspect = DataAspect())

    for i in 1:nc
        fi = G.cells.facePos[i]
        li = G.cells.facePos[i + 1] - 1
        cell_face_ids = G.cells.faces[fi:li]
        nf = length(cell_face_ids)

        # Collect ordered node indices for this cell polygon
        node_ids = Vector{Int}(undef, nf)
        for (j, fid) in enumerate(cell_face_ids)
            n1 = G.faces.nodes[G.faces.nodePos[fid]]
            n2 = G.faces.nodes[G.faces.nodePos[fid] + 1]
            node_ids[j] = G.faces.neighbors[fid, 1] == i ? n1 : n2
        end

        xs = coords[node_ids, 1]
        ys = coords[node_ids, 2]
        poly_points = Point2f.(xs, ys)

        if color_cells
            poly!(ax, poly_points;
                color = Float64(i),
                colormap = colormap,
                colorrange = (1.0, Float64(nc)),
                strokecolor = strokecolor,
                strokewidth = strokewidth,
            )
        else
            poly!(ax, poly_points;
                color = :transparent,
                strokecolor = strokecolor,
                strokewidth = strokewidth,
            )
        end
    end

    if show_nodes
        scatter!(ax, coords[:, 1], coords[:, 2];
            color = node_color,
            markersize = node_size,
        )
    end

    return fig
end

end # module PolytopeMeshesGLMakieExt
