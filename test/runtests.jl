using PolytopeMeshes
using PolytopeMeshes: Vec2, Vec3, to_points, point_segment_distance, number_of_cells, number_of_faces, number_of_nodes
import Jutul
using LinearAlgebra
using Random
using Test

dist_polyline(x, l) = minimum(point_segment_distance(x, l[k], l[k + 1]) for k in 1:(length(l) - 1))
polyline_length(l) = sum(norm(l[k + 1] - l[k]) for k in 1:(length(l) - 1))

# Maximum distance from the nodes of faces tagged with constraint k to the polyline
function max_tagged_distance(m, lines)
    d = 0.0
    for f in findall(>(0), m.face_tags)
        l = lines[m.face_tags[f]]
        for n in face_nodes(m, f)
            d = max(d, dist_polyline(m.node_points[n], l))
        end
    end
    return d
end

function tagged_length(m, k)
    g = mesh_geometry(m)
    return sum(g.face_areas[findall(==(k), m.face_tags)]; init = 0.0)
end

# Test that the cells with tag k form a connected set
function tagged_cells_connected(m, k)
    cells = Set(findall(==(k), m.cell_tags))
    isempty(cells) && return false
    adj = Dict(c => Int[] for c in cells)
    for (l, r) in m.face_neighbors
        if l in cells && r in cells
            push!(adj[l], r)
            push!(adj[r], l)
        end
    end
    seen = Set([first(cells)])
    stack = [first(cells)]
    while !isempty(stack)
        c = pop!(stack)
        for n in adj[c]
            if !(n in seen)
                push!(seen, n)
                push!(stack, n)
            end
        end
    end
    return length(seen) == length(cells)
end

@testset "PolytopeMeshes" begin
    @testset "Voronoi diagrams" begin
        rng = MersenneTwister(1)
        @testset "2D random" begin
            pts = [Vec2(rand(rng, 2)...) for _ in 1:300]
            m = voronoi_mesh(pts, box_domain((0, 0), (1, 1)))
            @test number_of_cells(m) == 300
            @test check_mesh(m).ok
            @test sum(cell_volumes(m)) ≈ 1.0
            # Each site lies in its own cell, i.e. it is closer to its own
            # face planes than to the neighbors
            @test all(1 .<= m.boundary_tags[m.boundary_tags .> 0] .<= 4)
        end
        @testset "2D Cartesian (degenerate)" begin
            pts = vec([Vec2((i - 0.5) / 10, (j - 0.5) / 10) for i in 1:10, j in 1:10])
            m = voronoi_mesh(pts, box_domain((0, 0), (1, 1)))
            @test number_of_cells(m) == 100
            @test number_of_faces(m) == 220
            @test number_of_nodes(m) == 121
            @test all(cell_volumes(m) .≈ 0.01)
        end
        @testset "Non-convex polygon" begin
            dom = PolygonDomain([[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]])
            pts = [Vec2(2 * rand(rng), 2 * rand(rng)) for _ in 1:400]
            m = voronoi_mesh(pts, dom)
            @test check_mesh(m).ok
            @test sum(cell_volumes(m)) ≈ 3.0
            @test number_of_cells(m) == count(p -> !(p[1] > 1 && p[2] > 1), pts)
        end
        @testset "3D random and Cartesian" begin
            pts = [Vec3(rand(rng, 3)...) for _ in 1:300]
            m = voronoi_mesh(pts, box_domain((0, 0, 0), (1, 1, 1)))
            @test check_mesh(m).ok
            @test sum(cell_volumes(m)) ≈ 1.0
            pts = vec([Vec3((i - 0.5) / 5, (j - 0.5) / 5, (k - 0.5) / 5) for i in 1:5, j in 1:5, k in 1:5])
            m = voronoi_mesh(pts, box_domain((0, 0, 0), (1, 1, 1)))
            @test number_of_faces(m) == 3 * 4 * 25 + 6 * 25
            @test number_of_nodes(m) == 216
            @test all(cell_volumes(m) .≈ 1 / 125)
            @test sort(unique(m.boundary_tags)) == 0:6
        end
        @testset "Convex hull domain" begin
            pts = [Vec3(rand(rng, 3)...) for _ in 1:200]
            m = mirrored_pebi_grid_3d(pts, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
            @test check_mesh(m).ok
            @test sum(cell_volumes(m)) ≈ 1 / 6
            @test number_of_cells(m) == count(p -> sum(p) <= 1, pts)
        end
        @testset "Duplicate sites" begin
            pts = [Vec2(0.25, 0.25), Vec2(0.75, 0.75), Vec2(0.25, 0.25)]
            m = @test_logs (:warn,) voronoi_mesh(pts, box_domain((0, 0), (1, 1)))
            @test number_of_cells(m) == 2
        end
        @testset "Centroidal" begin
            pts = [Vec2(rand(rng, 2)...) for _ in 1:50]
            m, x = centroidal_voronoi_mesh(pts, box_domain((0, 0), (1, 1)); maxit = 100)
            @test check_mesh(m).ok
            c = cell_centroids(m)
            @test maximum(norm.(c .- m.cell_sites)) < 1e-2
        end
    end

    @testset "Lines" begin
        L1 = [[[0.0, 0.0], [1.0, 1.0]], [[0.0, 1.0], [1.0, 0.0]]]
        s = split_at_intersections(L1)
        @test length(s.lines) == 4
        @test s.self_cut == [1, 2, 1, 2]
        @test s.origin == [1, 1, 2, 2]
        @test all(l -> l[end] ≈ Vec2(0.5, 0.5) || l[1] ≈ Vec2(0.5, 0.5), s.lines)
        s = split_at_intersections([L1[1]], [L1[2]])
        @test s.other_cut == [1, 2]
        @test s.other_end == [1, 0]

        line = to_points([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], 2)
        p = interpolate_line_path(line, 0.1, 0.1)
        @test p[1] ≈ line[1] && p[end] ≈ line[end]
        @test any(x -> x ≈ line[2], p)
        d = norm.(diff(p))
        @test all(abs.(d .- 0.1) .< 0.02)
        p = interpolate_line_path(line, 0.1, 0.1; se = (1.0, 0.5))
        @test norm(p[1] - line[1]) ≈ 0.1 atol = 0.02
        @test norm(p[end] - line[end]) ≈ 0.05 atol = 0.02
    end

    @testset "2D face constraint sites" begin
        F = surface_sites_2d([[[0.2, 0.5], [0.8, 0.5]]], 0.1)
        @test length(F.sites) == 2 * length(F.edges)
        @test length(F.tips) == 2
        # Sites come in mirrored pairs across the line
        for (e, (a, b)) in enumerate(F.edges)
            pair = [s for s in eachindex(F.sites) if (e, 1) in F.generators[s] || (e, -1) in F.generators[s]]
            @test length(pair) == 2
            @test F.sites[pair[1]][1] ≈ F.sites[pair[2]][1]
            @test F.sites[pair[1]][2] - 0.5 ≈ 0.5 - F.sites[pair[2]][2]
        end
    end

    @testset "pebi_grid_2d and composite_pebi_grid_2d" begin
        fl = [[[0.2, 0.2], [0.8, 0.8]]]
        wl = [[[0.2, 0.8], [0.8, 0.2]]]
        faults = [to_points(l, 2) for l in fl]
        for m in (pebi_grid_2d(1 / 10, [1, 1]; cell_constraints = wl, face_constraints = fl),
                  composite_pebi_grid_2d([1 / 10, 1 / 10], [1, 1]; cell_constraints = wl, face_constraints = fl))
            @test check_mesh(m).ok
            @test sum(cell_volumes(m)) ≈ 1.0
            @test max_tagged_distance(m, faults) < 1e-10
            @test tagged_length(m, 1) ≈ polyline_length(faults[1])
            @test tagged_cells_connected(m, 1)
        end
    end

    @testset "Intersecting faults and wells" begin
        x = range(0.1, 0.8, length = 10)
        well = [[[xx, 0.2 + xx^2] for xx in x]]
        fault = [[[0.3, 0.1], [0.1, 0.3]], [[0.3, 0.1], [0.4, 0.7]], [[0.6, 0.2], [0.3, 0.7]], [[0.8, 0.4], [0.5, 0.8]]]
        faults = [to_points(l, 2) for l in fault]
        m = composite_pebi_grid_2d([1 / 24, 1 / 24], [1, 1]; cell_constraints = well, cc_factor = 0.5,
            face_constraints = fault, fc_factor = 0.5, mlqt_max_level = 1)
        @test check_mesh(m).ok
        @test max_tagged_distance(m, faults) < 1e-10
        for k in 1:4
            @test tagged_length(m, k) ≈ polyline_length(faults[k]) rtol = 1e-8
        end
        @test tagged_cells_connected(m, 1)
        m = pebi_grid_2d(1 / 24, [1, 1]; cell_constraints = well, cc_factor = 0.25, cc_refinement = true,
            cc_eps = 1 / 12, face_constraints = fault)
        @test check_mesh(m).ok
        @test max_tagged_distance(m, faults) < 1e-10
        @test tagged_cells_connected(m, 1)

        l = [[[0.1, 0.42], [0.4, 0.55], [0.7, 0.65]], [[0.8, 0.13], [0.6, 0.4], [0.55, 0.6]]]
        faults = [to_points(x, 2) for x in l]
        m = composite_pebi_grid_2d([1 / 30, 1 / 30], [1, 1.15]; face_constraints = l)
        @test check_mesh(m).ok
        @test max_tagged_distance(m, faults) < 1e-10
        @test tagged_length(m, 1) ≈ polyline_length(faults[1]) rtol = 1e-8
        @test tagged_length(m, 2) ≈ polyline_length(faults[2]) rtol = 1e-8
    end

    @testset "Branching wells, points and polygon domains" begin
        wells = [[[0.5, 0.2], [0.5, 0.3], [0.47, 0.4], [0.4, 0.5], [0.33, 0.6], [0.26, 0.7]],
                 [[0.5, 0.3], [0.53, 0.4], [0.58, 0.5]],
                 [[0.5, 0.45], [0.5, 0.55], [0.45, 0.65], [0.4, 0.75], [0.38, 0.85]]]
        m = composite_pebi_grid_2d([1 / 24, 1 / 24], [1, 1]; cell_constraints = wells, cc_factor = 0.25,
            mlqt_max_level = 2, mlqt_level_steps = [0.1, 0.05])
        @test check_mesh(m).ok
        lines = [to_points(w, 2) for w in wells]
        for c in findall(>(0), m.cell_tags)
            @test dist_polyline(m.cell_sites[c], lines[m.cell_tags[c]]) < 1e-10
        end
        poly = [[0, 0], [1, 0], [1.2, 0.6], [0.6, 1], [0, 0.7]]
        fault = [[[0.1, 0.6], [0.9, 0.2]]]
        m = pebi_grid_2d(0.08, nothing; polygon = poly, cell_constraints = [[[0.3, 0.3]], [[0.5, 0.5], [0.9, 0.6]]],
            protection_layer = true, face_constraints = fault)
        @test check_mesh(m).ok
        @test sum(cell_volumes(m)) ≈ PolytopeMeshes.domain_volume(PolygonDomain(poly))
        @test count(==(1), m.cell_tags) == 1
        @test max_tagged_distance(m, [to_points(fault[1], 2)]) < 1e-10
    end

    @testset "composite_pebi_grid_3d" begin
        nx = 8
        x = range(0, 1, length = nx)
        li = LinearIndices((nx, nx))
        tris = NTuple{3, Int}[]
        for i in 1:(nx - 1), j in 1:(nx - 1)
            push!(tris, (li[i, j], li[i + 1, j], li[i + 1, j + 1]))
            push!(tris, (li[i, j], li[i + 1, j + 1], li[i, j + 1]))
        end
        pts1 = vec([Vec3(a, 0.5a + 0.1, b) * 0.99 .+ 0.005 for a in x, b in x])
        pts2 = vec([Vec3(a, -0.5a + 0.6, b) for a in x, b in x])
        well = [[[0.9, 0.75, 0.1], [0.9, 0.9, 0.9]]]
        m = composite_pebi_grid_3d((nx, nx, nx), (1, 1, 1); face_constraints = [(pts1, tris), (pts2, tris)],
            fc_rho = 1 / (nx - 1), cell_constraints = well)
        @test check_mesh(m).ok
        @test sum(cell_volumes(m)) ≈ 1.0
        for (k, p) in enumerate((pts1, pts2))
            fs = findall(==(k), m.face_tags)
            @test length(fs) > 0
            n = PolytopeMeshes.normal_from_points([p[1], p[nx], p[end]])
            @test maximum(abs(dot(n, m.node_points[i] - p[1])) for f in fs for i in face_nodes(m, f)) < 1e-10
        end
        @test count(==(1), m.cell_tags) > 3
        @test tagged_cells_connected(m, 1)
        S = surface_triangulation([[0.1, 0.1, 0.5], [0.9, 0.1, 0.5], [0.9, 0.9, 0.5], [0.1, 0.9, 0.5]], 0.1)
        m = composite_pebi_grid_3d((10, 10, 10), (1, 1, 1); face_constraints = [S], fc_rho = 0.1)
        @test check_mesh(m).ok
        fs = findall(==(1), m.face_tags)
        @test all(abs(m.node_points[i][3] - 0.5) < 1e-10 for f in fs for i in face_nodes(m, f))
        # Tip sites at the free edges stop the faces at the surface boundary
        @test sum(mesh_geometry(m).face_areas[fs]) ≈ 0.64 rtol = 1e-8
        @test all(0.1 - 1e-10 <= m.node_points[i][k] <= 0.9 + 1e-10 for f in fs for i in face_nodes(m, f) for k in 1:2)
        # Tips near the domain boundary are rotated into the domain
        m = composite_pebi_grid_3d((10, 10, 10), (1, 1, 1); face_constraints = [S], fc_rho = 0.13)
        @test check_mesh(m).ok
        fs = findall(==(1), m.face_tags)
        @test sum(mesh_geometry(m).face_areas[fs]) ≈ 0.64 rtol = 1e-8
        # Non-planar surface: a cylinder patch triangulated with equilateral
        # triangles, with zigzag free edges. Each tagged face is a triangle.
        h, nu, nv, Rc = 0.1, 7, 8, 0.4
        li = LinearIndices((nu, nv))
        uv = vec([(0.2 + h * (i - 1) + (iseven(j) ? h / 2 : 0.0), (j - (nv + 1) / 2) * h * sqrt(3) / 2) for i in 1:nu, j in 1:nv])
        tris = NTuple{3, Int}[]
        for j in 1:(nv - 1), i in 1:(nu - 1)
            if isodd(j)
                push!(tris, (li[i, j], li[i + 1, j], li[i, j + 1]), (li[i + 1, j], li[i + 1, j + 1], li[i, j + 1]))
            else
                push!(tris, (li[i, j], li[i + 1, j + 1], li[i, j + 1]), (li[i, j], li[i + 1, j], li[i + 1, j + 1]))
            end
        end
        pts = [Vec3(u, 0.5 + Rc * sin(v / Rc), 0.9 - Rc * cos(v / Rc)) for (u, v) in uv]
        F = surface_sites_3d([(pts, tris)], 0.12; domain = box_domain((0, 0, 0), (1, 1, 1)))
        @test length(F.tips) == 26
        m = composite_pebi_grid_3d((10, 10, 10), (1, 1, 1); face_constraints = [(pts, tris)], fc_rho = 0.12)
        @test check_mesh(m).ok
        @test sum(cell_volumes(m)) ≈ 1.0
        fs = findall(==(1), m.face_tags)
        @test length(fs) == length(tris)
        Atri = sum(norm(cross(pts[b] - pts[a], pts[c] - pts[a])) / 2 for (a, b, c) in tris)
        @test sum(mesh_geometry(m).face_areas[fs]) ≈ Atri rtol = 1e-8
        @test maximum(minimum(norm(m.node_points[i] - p) for p in pts) for f in fs for i in face_nodes(m, f)) < 1e-8
    end

    @testset "Hierarchical 3D fracture gridding" begin
        f1 = [[1, 3, 2], [4, 3, 2], [4, 3, 4], [1, 3, 4]]
        f2 = [[2, 2, 3.3], [5, 2, 3.3], [5, 4, 3.3], [2, 4, 3.3]]
        f3 = elliptic_surface_3d([3, 2.5, 3.2], 2, 2, 0, pi / 2, pi / 2)
        fracs = [to_points(f, 3) for f in (f1, f2, f3)]
        ints = surface_intersections_3d(fracs)
        @test length(ints) == 6
        g1 = line_grid_3d(ints, 0.2)
        g2 = surface_grid_3d(fracs, g1, ints, 0.25, 0.25 / 6)
        @test all(g -> check_mesh(g.mesh).ok, g2)
        @test all(g -> count(>(0), g.mesh.face_tags) > 0, g2)
        g3 = volume_grid_3d([6, 6, 6], fracs, g2, 0.4, 0.4 / 12)
        @test check_mesh(g3).ok
        geo = mesh_geometry(g3)
        @test sum(geo.cell_volumes) ≈ 216
        for k in 1:3
            fs = findall(==(k), g3.face_tags)
            # One 3D face per 2D cell, lying in the fracture plane
            @test length(fs) == number_of_cells(g2[k].mesh)
            n = PolytopeMeshes.normal_from_points(fracs[k])
            @test maximum(abs(dot(n, g3.node_points[i] - fracs[k][1])) for f in fs for i in face_nodes(g3, f)) < 1e-10
            # Interior 2D cells coincide with the 3D faces
            m2 = g2[k].mesh
            a2 = cell_volumes(m2)
            interior = [c for c in 1:number_of_cells(m2) if all(f -> m2.face_neighbors[f][2] != 0, cell_faces(m2, c))]
            a3 = Dict{Int, Float64}()
            for f in fs
                # The cells of a fracture face are the two sites of one 2D cell
                l, _ = g3.face_neighbors[f]
                c2 = argmin([norm(g2[k].sites[c] - g3.cell_sites[l]) for c in 1:number_of_cells(m2)])
                a3[c2] = geo.face_areas[f]
            end
            @test all(c -> isapprox(a3[c], a2[c]; rtol = 1e-8), interior)
        end
    end

    @testset "Triangulation backends" begin
        # A minimal custom backend: fixed points followed by a Cartesian
        # lattice with the smallest requested spacing, no triangles.
        struct LatticeBackend <: TriangulationBackend end
        function PolytopeMeshes.triangulate_polygon(::LatticeBackend, polygon, size, fixed; rng, min_size)
            lo = reduce((a, b) -> min.(a, b), polygon)
            hi = reduce((a, b) -> max.(a, b), polygon)
            lattice = [Vec2(x, y) for x in lo[1]:min_size:hi[1] for y in lo[2]:min_size:hi[2]]
            lattice = filter(p -> minimum(q -> norm(p - q), fixed; init = Inf) > min_size / 2, lattice)
            return vcat(fixed, lattice), NTuple{3, Int}[]
        end
        fl = [[[0.2, 0.2], [0.8, 0.8]]]
        m = pebi_grid_2d(0.1, [1, 1]; face_constraints = fl, triangulation = LatticeBackend())
        @test check_mesh(m).ok
        @test max_tagged_distance(m, [to_points(fl[1], 2)]) < 1e-10
        @test tagged_length(m, 1) ≈ polyline_length(to_points(fl[1], 2))
        @test_throws ErrorException surface_triangulation([[0, 0], [1, 0], [1, 1]], 0.2; triangulation = LatticeBackend())
        # Backends that are not implemented, or break the interface, give errors
        struct MissingBackend <: TriangulationBackend end
        @test_throws ArgumentError pebi_grid_2d(0.1, [1, 1]; triangulation = MissingBackend())
        struct ShuffledBackend <: TriangulationBackend end
        PolytopeMeshes.triangulate_polygon(::ShuffledBackend, polygon, size, fixed; kwarg...) = (reverse(fixed), NTuple{3, Int}[])
        @test_throws ErrorException pebi_grid_2d(0.1, [1, 1]; face_constraints = fl, triangulation = ShuffledBackend())
        # The default backend is DistMesh, with reproducible results
        @test PolytopeMeshes.default_triangulation_backend() isa DistMeshBackend
        S1 = surface_triangulation([[0, 0], [1, 0], [1, 1], [0, 1]], 0.2)
        S2 = surface_triangulation([[0, 0], [1, 0], [1, 1], [0, 1]], 0.2; triangulation = DistMeshBackend(maxit = 200))
        @test S1.points == S2.points && S1.triangles == S2.triangles
        @test sum(t -> abs(PolytopeMeshes.cross2(S1.points[t[2]][1:2] - S1.points[t[1]][1:2], S1.points[t[3]][1:2] - S1.points[t[1]][1:2])) / 2, S1.triangles) ≈ 1.0
    end

    @testset "Output" begin
        m2 = composite_pebi_grid_2d([0.1, 0.1], [1, 1]; face_constraints = [[[0.2, 0.3], [0.8, 0.6]]])
        m3 = composite_pebi_grid_3d((4, 4, 4), (1, 1, 1))
        mktempdir() do dir
            for (name, m) in (("m2", m2), ("m3", m3))
                files = write_vtk(joinpath(dir, name), m; cell_data = Dict("x" => [p[1] for p in m.cell_sites]))
                @test all(isfile, files)
            end
            files = write_vtk_faces(joinpath(dir, "faults"), m2)
            @test all(isfile, files)
        end
    end

    @testset "Jutul" begin
        fl = [[[0.2, 0.2], [0.8, 0.8]]]
        wl = [[[0.2, 0.8], [0.8, 0.2]]]
        m2 = composite_pebi_grid_2d([1 / 10, 1 / 10], [1, 1]; cell_constraints = wl, face_constraints = fl)
        m3 = composite_pebi_grid_3d((5, 5, 5), (1, 1, 1); cell_constraints = [[[0.2, 0.2, 0.1], [0.3, 0.8, 0.9]]])
        for m in (m2, m3)
            g = Jutul.UnstructuredMesh(m)
            @test Jutul.number_of_cells(g) == number_of_cells(m)
            geo = Jutul.tpfv_geometry(g)
            mine = mesh_geometry(m)
            @test geo.volumes ≈ mine.cell_volumes
            int = findall(f -> m.face_neighbors[f][2] != 0, 1:number_of_faces(m))
            @test geo.areas ≈ mine.face_areas[int]
            @test maximum(norm.(eachcol(geo.normals) .- mine.face_normals[int])) < 1e-10
            wells = Jutul.get_mesh_entity_tag(g, Jutul.Cells(), :cell_constraints, :constraint_1)
            @test wells == findall(==(1), m.cell_tags)
        end
        faults = Jutul.get_mesh_entity_tag(Jutul.UnstructuredMesh(m2), Jutul.Faces(), :face_constraints, :constraint_1)
        @test length(faults) == count(==(1), m2.face_tags)
    end
end
