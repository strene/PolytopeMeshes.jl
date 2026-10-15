# Static mesh plots.

_point(p::SVector{2}) = Point2f(p[1], p[2])
_point(p::SVector{3}) = Point3f(p[1], p[2], p[3])

function _triangle_matrix(tris)
    T = zeros(Int, length(tris), 3)
    for (i, t) in enumerate(tris)
        T[i, 1], T[i, 2], T[i, 3] = t
    end
    return T
end

_segments(edges) = [_point(p) for e in edges for p in e]

function _new_axis(pos, ::PolytopalMesh{2}; kwarg...)
    return Axis(pos; aspect = DataAspect(), kwarg...)
end

function _new_axis(pos, ::PolytopalMesh{3}; kwarg...)
    return Axis3(pos; aspect = :data, kwarg...)
end

# One distinct color per constraint tag
function _tag_colors(tags, colormap)
    palette = Makie.to_colormap(colormap)
    return [palette[mod1(t, length(palette))] for t in tags]
end

function _cell_selection(m::PolytopalMesh, cells)
    nc = number_of_cells(m)
    isnothing(cells) && return trues(nc)
    sel = falses(nc)
    if eltype(cells) == Bool
        sel .= cells
    else
        sel[collect(cells)] .= true
    end
    return sel
end

function PolytopeMeshes.plot_mesh_impl(m::PolytopalMesh; size = (1000, 800), colorbar = true, axis = (;), kwarg...)
    fig = Figure(; size = size)
    ax = _new_axis(fig[1, 1], m; axis...)
    p = PolytopeMeshes.plot_mesh!(ax, m; kwarg...)
    color = get(kwarg, :color, nothing)
    if colorbar && color isa AbstractVector{<:Real}
        Colorbar(fig[1, 2], p)
    end
    return Makie.FigureAxisPlot(fig, ax, p)
end

function PolytopeMeshes.plot_mesh_impl!(ax, m::PolytopalMesh{D};
        cells = nothing,
        color = :lightblue,
        colormap = :viridis,
        colorrange = Makie.automatic,
        edges = true,
        edgecolor = :black,
        edgewidth = 0.5,
        constraints = true,
        constraint_colormap = :tab10,
        kwarg...
    ) where {D}
    sel = _cell_selection(m, cells)
    tri = plot_triangulation(m, sel)
    pts = _point.(tri.points)
    faces = _triangle_matrix(tri.triangles)
    if color isa AbstractVector
        length(color) == number_of_cells(m) || throw(ArgumentError("Expected one color value per cell"))
        vcolor = color[tri.vertex_cell]
        if eltype(vcolor) <: Real
            vcolor = Float32.(vcolor)
        end
    else
        vcolor = color
    end
    shading = D == 2 ? NoShading : Makie.automatic
    p = mesh!(ax, pts, faces; color = vcolor, colormap = colormap, colorrange = colorrange, shading = shading, kwarg...)
    if edges && !isempty(tri.edges)
        linesegments!(ax, _segments(tri.edges); color = edgecolor, linewidth = edgewidth)
    end
    if constraints
        _plot_constraints!(ax, m, sel, constraint_colormap)
    end
    return p
end

function _plot_constraints!(ax, m::PolytopalMesh, sel, colormap)
    faces = [f for f in findall(>(0), m.face_tags) if any(c -> c > 0 && sel[c], m.face_neighbors[f])]
    pf = isempty(faces) ? nothing : PolytopeMeshes.plot_faces!(ax, m; faces = faces, colormap = colormap)
    cells = [c for c in findall(>(0), m.cell_tags) if sel[c]]
    pc = nothing
    if !isempty(cells)
        pc = scatter!(ax, _point.(m.cell_sites[cells]); color = _tag_colors(m.cell_tags[cells], colormap),
            markersize = 8, strokecolor = :black, strokewidth = 0.5)
    end
    return pf, pc
end

function PolytopeMeshes.plot_faces_impl!(ax, m::PolytopalMesh{D};
        faces = findall(>(0), m.face_tags),
        color = nothing,
        colormap = :tab10,
        linewidth = D == 2 ? 3 : 1,
        visible = true,
        kwarg...
    ) where {D}
    faces = eltype(faces) == Bool ? findall(faces) : collect(faces)
    isempty(faces) && return nothing
    fcolor = isnothing(color) ? _tag_colors(max.(m.face_tags[faces], 1), colormap) : color
    if D == 2
        pts = Point2f[]
        for f in faces
            a, b = face_nodes(m, f)
            push!(pts, _point(m.node_points[a]), _point(m.node_points[b]))
        end
        c = fcolor isa AbstractVector ? repeat(fcolor; inner = 2) : fcolor
        return linesegments!(ax, pts; color = c, linewidth = linewidth, visible = visible, kwarg...)
    else
        pts = Point3f[]
        tris = NTuple{3, Int}[]
        vcolor = fcolor isa AbstractVector ? similar(fcolor, 0) : fcolor
        seg = Point3f[]
        for (i, f) in enumerate(faces)
            nodes = face_nodes(m, f)
            n0 = length(pts)
            k = length(nodes)
            for n in nodes
                push!(pts, _point(m.node_points[n]))
                fcolor isa AbstractVector && push!(vcolor, fcolor[i])
            end
            for j in 2:(k - 1)
                push!(tris, (n0 + 1, n0 + j, n0 + j + 1))
            end
            for j in 1:k
                push!(seg, pts[n0 + j], pts[n0 + mod1(j + 1, k)])
            end
        end
        p = mesh!(ax, pts, _triangle_matrix(tris); color = vcolor, visible = visible, kwarg...)
        linesegments!(ax, seg; color = :black, linewidth = linewidth, visible = visible)
        return p
    end
end
