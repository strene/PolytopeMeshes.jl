# Interactive mesh viewer. Works with any Makie backend, but is meant for
# GLMakie, where the widgets are live.

function _cell_fields(m::PolytopalMesh, data)
    nc = number_of_cells(m)
    geo = mesh_geometry(m)
    fields = Pair{String, Vector{Float64}}[]
    push!(fields, "Cell constraints" => Float64.(m.cell_tags))
    push!(fields, (length(first(m.node_points)) == 2 ? "Area" : "Volume") => geo.cell_volumes)
    push!(fields, "Cell index" => collect(Float64, 1:nc))
    for (k, v) in pairs(data)
        length(v) == nc || throw(ArgumentError("Field $k has $(length(v)) values, expected one per cell ($nc)"))
        push!(fields, string(k) => Float64.(v))
    end
    return fields
end

function _range(v)
    lo, hi = extrema(v)
    if !(hi > lo)
        lo, hi = lo - 0.5, hi + 0.5
    end
    return (Float32(lo), Float32(hi))
end

function PolytopeMeshes.plot_polymesh_interactive_impl(m::PolytopalMesh{D};
        data = Dict{String, Vector{Float64}}(),
        size = (1400, 900),
        colormap = :viridis,
        constraint_colormap = :tab10
    ) where {D}
    nc = number_of_cells(m)
    fields = _cell_fields(m, data)
    names = first.(fields)
    fig = Figure(; size = size)
    ax = _new_axis(fig[1, 1], m)
    panel = GridLayout(fig[1, 3]; tellheight = false, valign = :top)
    row = Ref(0)
    nextrow() = (row[] += 1)

    Label(panel[nextrow(), 1:2], "$(nc) cells, $(number_of_faces(m)) faces"; font = :bold, halign = :left)
    Label(panel[nextrow(), 1:2], "Color by"; halign = :left)
    menu = Menu(panel[nextrow(), 1:2]; options = names, default = names[min(2, length(names))], width = 270)

    toggle_row(label, active) = begin
        r = nextrow()
        t = Toggle(panel[r, 1]; active = active)
        Label(panel[r, 2], label; halign = :left)
        t
    end
    t_cells = toggle_row("Cells", true)
    t_edges = toggle_row("Cell edges", true)
    t_faces = toggle_row("Face constraints", true)
    t_sites = toggle_row("Cell constraint sites", true)
    t_only = toggle_row("Only cell constraint cells", false)

    # Cell selection from the site positions
    sites = m.cell_sites
    sliders = []
    if D == 3
        for (i, axname) in enumerate(("x", "y", "z"))
            lo, hi = extrema(s[i] for s in sites)
            hi > lo || (hi = lo + 1)
            r = nextrow()
            sl = IntervalSlider(panel[r + 1, 1:2]; range = range(lo, hi, length = 201), startvalues = (lo, hi), width = 270)
            txt = lift(iv -> "Sites in $axname: " * string(round(iv[1], sigdigits = 3)) * " – " * string(round(iv[2], sigdigits = 3)), sl.interval)
            Label(panel[r, 1:2], txt; halign = :left)
            nextrow()
            push!(sliders, sl)
        end
    end
    selection = if D == 3
        lift(sliders[1].interval, sliders[2].interval, sliders[3].interval, t_only.active) do ix, iy, iz, only
            ints = (ix, iy, iz)
            sel = [all(ints[i][1] <= s[i] <= ints[i][2] for i in 1:3) for s in sites]
            only && (sel .&= m.cell_tags .> 0)
            sel
        end
    else
        lift(only -> only ? m.cell_tags .> 0 : trues(nc), t_only.active)
    end

    # Cells
    tri = lift(sel -> plot_triangulation(m, sel), selection)
    field = lift(name -> last(fields[findfirst(==(name), names)]), menu.selection)
    # Plot buffers for a selection. An empty selection gives a degenerate
    # triangle that is not drawn.
    function cell_buffers(t, v)
        if isempty(t.triangles)
            p0 = _point(m.node_points[1])
            return [p0, p0, p0], [1 2 3], Float32[v[1], v[1], v[1]]
        end
        return _point.(t.points), _triangle_matrix(t.triangles), Float32.(v[t.vertex_cell])
    end
    p0, f0, c0 = cell_buffers(tri[], field[])
    edge_pts = Observable(_segments(tri[].edges))
    shading = D == 2 ? NoShading : Makie.automatic
    pmesh = mesh!(ax, p0, f0; color = c0, colormap = colormap, colorrange = _range(field[]),
        shading = shading, visible = t_cells.active)
    onany(tri, field) do t, v
        # Update positions, faces and colors together, since their sizes change
        newp, newf, newc = cell_buffers(t, v)
        Makie.update!(pmesh, newp, newf; color = newc, colorrange = _range(v))
        edge_pts[] = isempty(t.edges) ? Makie.Point{D, Float32}[] : _segments(t.edges)
    end
    linesegments!(ax, edge_pts; color = :black, linewidth = 0.5,
        visible = lift((a, b) -> a && b, t_edges.active, t_cells.active))
    Colorbar(fig[1, 2], pmesh)

    # Constraints (all of them, independent of the selection)
    if any(>(0), m.face_tags)
        PolytopeMeshes.plot_polymesh_faces!(ax, m; colormap = constraint_colormap, visible = t_faces.active)
    end
    wells = findall(>(0), m.cell_tags)
    if !isempty(wells)
        sc = scatter!(ax, _point.(m.cell_sites[wells]); color = _tag_colors(m.cell_tags[wells], constraint_colormap),
            markersize = 8, strokecolor = :black, strokewidth = 0.5, visible = t_sites.active)
    end

    # Camera presets
    if D == 3
        r = nextrow()
        views = [("3D", 1.275pi, pi / 8), ("xy", -pi / 2, pi / 2), ("xz", -pi / 2, 0.0), ("yz", 0.0, 0.0)]
        bl = GridLayout(panel[r, 1:2])
        for (i, (name, az, el)) in enumerate(views)
            b = Button(bl[1, i]; label = name)
            on(b.clicks) do _
                ax.azimuth[] = az
                ax.elevation[] = el
            end
        end
    end
    r = nextrow()
    reset = Button(panel[r, 1:2]; label = "Reset view")
    on(reset.clicks) do _
        reset_limits!(ax)
    end
    colsize!(fig.layout, 3, Fixed(300))
    PolytopeMeshes.display_interactive(fig, Val(PolytopeMeshes.current_backend_name()))
    return fig
end
