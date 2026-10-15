# Plotting interface. The functions below are implemented in package
# extensions: PolytopeMeshesMakieExt (any Makie backend, e.g. CairoMakie for
# static figures) and PolytopeMeshesGLMakieExt (interactive windows). This
# follows the approach of Jutul and JutulDarcy.

"""
    check_plotting_availability(; throw = true, interactive = false)

Check that a Makie backend is loaded (e.g. `using GLMakie` or `using
CairoMakie`). With `throw = false` the availability is returned as a boolean.
If `interactive = true`, warn if the backend is not interactive (GLMakie is
recommended for [`plot_interactive`](@ref)).
"""
function check_plotting_availability(; throw = true, interactive = false)
    ok = try
        check_plotting_availability_impl()
    catch e
        if e isa MethodError
            throw && error("Plotting is not available. Load a Makie backend first, e.g. `using GLMakie` (interactive, recommended) or `using CairoMakie` (static figures).")
            false
        else
            rethrow(e)
        end
    end
    if ok && interactive && !plotting_is_interactive()
        @warn "Interactive plotting works best with GLMakie (`using GLMakie`). The current backend gives a static figure."
    end
    return ok
end

function check_plotting_availability_impl end

# Name of the active Makie backend (implemented in the Makie extension)
function current_backend_name end

# Interactive backends (GLMakie) add methods for their name
backend_is_interactive(::Val) = false
plotting_is_interactive() = backend_is_interactive(Val(current_backend_name()))

# Show a figure from plot_interactive. GLMakie opens a new window.
display_interactive(fig, ::Val) = display(fig)

"""
    plot_mesh(mesh; kwarg...)

Plot a [`PolytopalMesh`](@ref) in a new figure and return `(fig, ax, plt)`.
The cells are colored by `color` (a color, or one value per cell), and the
constraints are highlighted: faces on face constraints (faults/fractures)
and cells on cell constraints (wells). Requires a Makie backend.

# Keyword arguments
- `cells = nothing`: indices or boolean mask of the cells to plot (all by
  default). In 3D only the outer faces of the selected cells are drawn.
- `color = :lightblue`: color of the cells, or a vector of cell values.
- `colormap = :viridis`, `colorrange`: used when `color` is a vector.
- `edges = true`, `edgecolor = :black`, `edgewidth = 0.5`: cell edges.
- `constraints = true`: highlight constraint faces and cells.
- `constraint_colormap = :tab10`: colors of the constraints by tag.
- `colorbar = true`: add a colorbar when `color` is a vector.
- `size = (1000, 800)`: figure size.
"""
function plot_mesh(arg...; kwarg...)
    check_plotting_availability()
    return plot_mesh_impl(arg...; kwarg...)
end

"""
    plot_mesh!(ax, mesh; kwarg...)

Plot a mesh into an existing Makie axis (`Axis` in 2D, `Axis3` or `LScene` in
3D). Takes the same keyword arguments as [`plot_mesh`](@ref), except
`colorbar` and `size`, and returns the plot object of the cells.
"""
function plot_mesh!(arg...; kwarg...)
    check_plotting_availability()
    return plot_mesh_impl!(arg...; kwarg...)
end

"""
    plot_cell_data(mesh, data; kwarg...)

Plot one value per cell on the mesh. Equivalent to
`plot_mesh(mesh; color = data, kwarg...)`.
"""
plot_cell_data(mesh, data; kwarg...) = plot_mesh(mesh; color = data, kwarg...)

"""
    plot_cell_data!(ax, mesh, data; kwarg...)

Plot one value per cell on the mesh into an existing axis.
"""
plot_cell_data!(ax, mesh, data; kwarg...) = plot_mesh!(ax, mesh; color = data, kwarg...)

"""
    plot_faces!(ax, mesh; faces = findall(>(0), mesh.face_tags), kwarg...)

Plot a subset of the faces (by default those on face constraints) into an
existing axis, as lines in 2D and polygons in 3D, colored by `color` (a
color or one value per face in `faces`).
"""
function plot_faces!(arg...; kwarg...)
    check_plotting_availability()
    return plot_faces_impl!(arg...; kwarg...)
end

"""
    plot_interactive(mesh; data = Dict(), kwarg...)

Interactive viewer for a mesh (GLMakie recommended). The cells can be colored
by the cell tags, cell volumes, cell index or any of the cell fields given
in `data` (a dictionary or named tuple of per-cell vectors), selected from a
menu. Toggles show/hide cell edges, constraint faces and the domain.

In 3D, interval sliders select the cells by the position of their sites in
`x`, `y` and `z`, which makes it possible to look inside the mesh. A toggle
restricts the view to cells on cell constraints, and buttons set the camera
to the coordinate planes.

Returns the figure. With GLMakie the figure opens in its own window.
"""
function plot_interactive(arg...; kwarg...)
    check_plotting_availability(interactive = true)
    return plot_interactive_impl(arg...; kwarg...)
end

function plot_mesh_impl end
function plot_mesh_impl! end
function plot_faces_impl! end
function plot_interactive_impl end

# Makie-independent geometry for plotting ---------------------------------

"""
    plot_triangulation(mesh, cells = 1:number_of_cells(mesh))

Triangulate the selected cells for plotting. In 2D the cells are fanned from
their first node; in 3D the faces on the outside of the selection are fanned.
The vertices are duplicated for each cell so that cell values give flat
colors.

Returns a named tuple with `points`, `triangles` (`NTuple{3, Int}`),
`vertex_cell` (cell of each vertex) and `edges` (pairs of points for the
cell edges that are visible).
"""
function plot_triangulation(m::PolytopalMesh{D}, cells = 1:number_of_cells(m)) where {D}
    nc = number_of_cells(m)
    selected = falses(nc)
    if eltype(cells) == Bool
        length(cells) == nc || throw(ArgumentError("Boolean cell mask must have one entry per cell"))
        selected .= cells
    else
        selected[collect(cells)] .= true
    end
    T = SVector{D, Float64}
    points = T[]
    triangles = NTuple{3, Int}[]
    vertex_cell = Int[]
    edges = Tuple{T, T}[]
    function add_loop!(loop, c)
        n0 = length(points)
        k = length(loop)
        for n in loop
            push!(points, m.node_points[n])
            push!(vertex_cell, c)
        end
        for i in 2:(k - 1)
            push!(triangles, (n0 + 1, n0 + i, n0 + i + 1))
        end
    end
    if D == 2
        for c in 1:nc
            selected[c] || continue
            add_loop!(cell_nodes(m, c), c)
        end
        for f in 1:number_of_faces(m)
            l, r = m.face_neighbors[f]
            if selected[l] || (r > 0 && selected[r])
                a, b = face_nodes(m, f)
                push!(edges, (m.node_points[a], m.node_points[b]))
            end
        end
    else
        for f in 1:number_of_faces(m)
            l, r = m.face_neighbors[f]
            sl = selected[l]
            sr = r > 0 && selected[r]
            sl == sr && continue
            nodes = face_nodes(m, f)
            # Orient the face outwards from the selected cell
            loop = sl ? collect(nodes) : reverse(nodes)
            add_loop!(loop, sl ? l : r)
            k = length(loop)
            for i in 1:k
                push!(edges, (m.node_points[loop[i]], m.node_points[loop[mod1(i + 1, k)]]))
            end
        end
    end
    return (points = points, triangles = triangles, vertex_cell = vertex_cell, edges = edges)
end
