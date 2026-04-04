using Test
using PolytopeMeshes

@testset "PolytopeMeshes.jl" begin

    @testset "Mesh data structure" begin
        # Create a simple 2-cell mesh (two triangles)
        coords = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.5]
        cell_faces = [1, 2, 3, 4, 5, 3]
        cell_facePos = [1, 4, 7]
        face_nodes = [1, 2, 2, 5, 5, 1, 2, 3, 3, 5]
        face_nodePos = [1, 3, 5, 7, 9, 11]
        face_neighbors = [1 0; 1 2; 1 0; 2 0; 2 0]

        G = UnstructuredMesh(coords, cell_faces, cell_facePos, face_nodes,
            face_nodePos, face_neighbors; meshdim=2)

        @test G.cells.num == 2
        @test G.faces.num == 5
        @test G.nodes.num == 5
        @test G.meshdim == 2
    end

    @testset "Point in polygon" begin
        poly = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]

        @test PolytopeMeshes.point_in_polygon(0.5, 0.5, poly) == true
        @test PolytopeMeshes.point_in_polygon(1.5, 0.5, poly) == false
        @test PolytopeMeshes.point_in_polygon(-0.1, 0.5, poly) == false
    end

    @testset "Line-line intersection" begin
        L1 = [0.0 0.0 1.0 1.0]  # diagonal line
        L2 = [0.0 1.0 1.0 0.0]  # opposite diagonal

        X, Y, seg = PolytopeMeshes.line_line_intersection(L1, L2)
        @test seg[1, 1] == true
        @test isapprox(X[1, 1], 0.5; atol=1e-10)
        @test isapprox(Y[1, 1], 0.5; atol=1e-10)

        # Parallel lines - no intersection
        L1 = [0.0 0.0 1.0 0.0]
        L2 = [0.0 1.0 1.0 1.0]
        _, _, seg2 = PolytopeMeshes.line_line_intersection(L1, L2)
        @test seg2[1, 1] == false
    end

    @testset "Remove conflict points" begin
        P1 = [1.0 1.0; 2.0 2.0; 5.0 5.0; 8.0 8.0]
        P2 = [1.0 1.0; 5.0 5.0]
        dist = [1.5, 1.5]

        result, removed = remove_conflict_points(P1, P2, dist)
        @test removed[1] == true   # (1,1) is within 1.5 of P2[1]
        @test removed[2] == true   # (2,2) is within 1.5 of P2[1]
        @test removed[3] == true   # (5,5) is within 1.5 of P2[2]
        @test removed[4] == false  # (8,8) is far from both
    end

    @testset "Clip polygon" begin
        poly = [0.0 0.0; 2.0 0.0; 2.0 2.0; 0.0 2.0]
        normals = [1.0 0.0]  # clip at x = 1
        points = [1.0 1.0]

        result = PolytopeMeshes.clip_polygon(poly, normals, points)
        @test size(result, 1) == 4  # Should be a rectangle
        @test all(result[:, 1] .<= 1.0 + 1e-10)  # All x <= 1
    end

    @testset "Interpolation" begin
        path = [0.0 0.0; 1.0 0.0; 1.0 1.0]
        p, dt = PolytopeMeshes.eq_interpret(path, 0.5, [0.0, 0.0])
        @test size(p, 1) >= 3
        @test isapprox(p[1, :], [0.0, 0.0]; atol=1e-10)
    end

    @testset "MLQT refinement" begin
        center = [0.5, 0.5]
        bndr = [0.5 0.5]
        cell_size = [1.0, 1.0]

        result = PolytopeMeshes.mlqt(center, bndr, cell_size;
            level=1, max_level=1)
        @test length(result) == 4  # Refined into 4 quadrants

        # No refinement when far away
        center2 = [10.0, 10.0]
        result2 = PolytopeMeshes.mlqt(center2, bndr, cell_size;
            level=1, max_level=1)
        @test length(result2) == 1  # Not refined
    end

    @testset "Delaunay triangulation" begin
        # Simple square of points
        pts = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        tri = PolytopeMeshes._delaunay_2d(pts)
        @test size(tri, 1) == 2  # Square -> 2 triangles
        @test size(tri, 2) == 3
    end

    @testset "Triangulation backend" begin
        # Default backend should be BowyerWatson
        @test default_triangulation_backend() isa BowyerWatson

        # triangulate_2d with BowyerWatson matches _delaunay_2d
        pts = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        tri_direct = PolytopeMeshes._delaunay_2d(pts)
        tri_backend = triangulate_2d(BowyerWatson(), pts)
        @test tri_direct == tri_backend

        # Edge case: fewer than 3 points
        @test size(triangulate_2d(BowyerWatson(), zeros(2, 2)), 1) == 0

        # set/get default backend round-trips
        original = default_triangulation_backend()
        new_bw = BowyerWatson()
        set_default_triangulation_backend!(new_bw)
        @test default_triangulation_backend() === new_bw
        set_default_triangulation_backend!(original)

        # clipped_voronoi_2d accepts triangulation kwarg
        sites = zeros(9, 2)
        idx = 1
        for y in [0.25, 0.5, 0.75]
            for x in [0.25, 0.5, 0.75]
                sites[idx, :] = [x, y]
                idx += 1
            end
        end
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        G = clipped_voronoi_2d(sites, bnd; triangulation=BowyerWatson())
        @test G.cells.num == 9
        @test G.meshdim == 2
    end

    @testset "Clipped Voronoi 2D - simple" begin
        # Regular mesh of points
        pts = zeros(9, 2)
        idx = 1
        for y in [0.25, 0.5, 0.75]
            for x in [0.25, 0.5, 0.75]
                pts[idx, :] = [x, y]
                idx += 1
            end
        end
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]

        G = clipped_voronoi_2d(pts, bnd)

        @test G.cells.num == 9
        @test G.faces.num > 0
        @test G.nodes.num > 0
        @test G.meshdim == 2

        # Compute geometry
        compute_geometry!(G)
        @test length(G.cells.volumes) == 9
        @test all(G.cells.volumes .> 0)
        @test isapprox(sum(G.cells.volumes), 1.0; atol=0.05)
    end

    @testset "Line sites 2D" begin
        cc = [Float64[0.2 0.2; 0.8 0.8]]
        gs = 0.1

        cc_pts, c_gs, _, _ = line_sites_2d(cc, gs)
        @test size(cc_pts, 1) > 0
        @test size(cc_pts, 2) == 2
        @test length(c_gs) == size(cc_pts, 1)
    end

    @testset "Surface sites 2D" begin
        fc = [Float64[0.2 0.5; 0.8 0.5]]
        gs = 0.2

        F = surface_sites_2d(fc, gs)
        @test size(F.f_pts, 1) > 0
        @test size(F.f_pts, 2) == 2
        @test size(F.c_CC, 1) > 0
        @test length(F.c_R) > 0
    end

    @testset "Surface sufficient condition" begin
        F = surface_sites_2d([Float64[0.2 0.5; 0.8 0.5]], 0.2)
        pts = rand(200, 2)

        filtered, removed = surface_suf_cond_2d(pts, F)
        @test size(filtered, 1) + count(removed) == 200
    end

    @testset "Composite Voronoi mesh 2D - basic" begin
        G, pts, F = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0])

        @test G.cells.num > 0
        @test G.faces.num > 0
        @test G.nodes.num > 0

        compute_geometry!(G)
        @test isapprox(sum(G.cells.volumes), 1.0; atol=0.1)
    end

    @testset "Composite Voronoi mesh 2D - with cell constraints" begin
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
            cell_constraints=wl)

        @test G.cells.num > 0
        @test any(G.cells.tag)  # Some cells should be tagged as wells
    end

    @testset "Composite Voronoi mesh 2D - with face constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        G, pts, F = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
            face_constraints=fl)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0
    end

    @testset "Composite Voronoi mesh 2D - with both constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
            cell_constraints=wl,
            face_constraints=fl)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0

        compute_geometry!(G)
        @test all(G.cells.volumes .> 0)
    end

    @testset "Extrude - single layer" begin
        # Create a simple 2D mesh and extrude it
        pts = zeros(9, 2)
        idx = 1
        for y in [0.25, 0.5, 0.75]
            for x in [0.25, 0.5, 0.75]
                pts[idx, :] = [x, y]
                idx += 1
            end
        end
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        G2d = clipped_voronoi_2d(pts, bnd)
        compute_geometry!(G2d)

        G3d = extrude(G2d, 1.0, 1)

        @test G3d.meshdim == 3
        @test G3d.nodes.num == G2d.nodes.num * 2
        @test G3d.cells.num == G2d.cells.num
        @test size(G3d.nodes.coords, 2) == 3

        # Compute 3D geometry
        compute_geometry!(G3d)
        @test all(G3d.cells.volumes .> 0)
        # Volume should equal 2D area * height (height = 1.0)
        @test isapprox(sum(G3d.cells.volumes), sum(G2d.cells.volumes) * 1.0; atol=0.05)
    end

    @testset "Extrude - multiple layers" begin
        pts = zeros(4, 2)
        pts[1, :] = [0.25, 0.25]
        pts[2, :] = [0.75, 0.25]
        pts[3, :] = [0.75, 0.75]
        pts[4, :] = [0.25, 0.75]
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        G2d = clipped_voronoi_2d(pts, bnd)
        compute_geometry!(G2d)

        # 3 layers from 0 to 1.0
        G3d = extrude(G2d, 1.0, 3)

        @test G3d.meshdim == 3
        @test G3d.cells.num == G2d.cells.num * 3
        @test G3d.nodes.num == G2d.nodes.num * 4  # 4 levels

        compute_geometry!(G3d)
        @test all(G3d.cells.volumes .> 0)
        @test isapprox(sum(G3d.cells.volumes), sum(G2d.cells.volumes) * 1.0; atol=0.05)
    end

    @testset "Extrude - variable layers" begin
        pts = zeros(4, 2)
        pts[1, :] = [0.25, 0.25]
        pts[2, :] = [0.75, 0.25]
        pts[3, :] = [0.75, 0.75]
        pts[4, :] = [0.25, 0.75]
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        G2d = clipped_voronoi_2d(pts, bnd)
        compute_geometry!(G2d)

        # 2 layers from 0 to 1, then 3 layers from 1 to 3
        G3d = extrude(G2d, [1.0, 3.0], [2, 3])

        @test G3d.meshdim == 3
        @test G3d.cells.num == G2d.cells.num * 5
        @test G3d.nodes.num == G2d.nodes.num * 6  # 6 levels

        compute_geometry!(G3d)
        @test all(G3d.cells.volumes .> 0)
        @test isapprox(sum(G3d.cells.volumes), sum(G2d.cells.volumes) * 3.0; atol=0.05)
    end

    @testset "Extrude - dim kwarg" begin
        pts = zeros(4, 2)
        pts[1, :] = [0.25, 0.25]
        pts[2, :] = [0.75, 0.25]
        pts[3, :] = [0.75, 0.75]
        pts[4, :] = [0.25, 0.75]
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        G2d = clipped_voronoi_2d(pts, bnd)

        # Extrude in x direction
        Gx = extrude(G2d, 2.0, 1; dim=:x)
        @test Gx.meshdim == 3
        # x-coords should contain the extrusion values (0 and 2)
        @test minimum(Gx.nodes.coords[:, 1]) ≈ 0.0
        @test maximum(Gx.nodes.coords[:, 1]) ≈ 2.0
        # y-coords should contain original x values
        @test minimum(Gx.nodes.coords[:, 2]) ≈ 0.0
        @test maximum(Gx.nodes.coords[:, 2]) ≈ 1.0

        # Extrude in y direction
        Gy = extrude(G2d, 2.0, 1; dim=:y)
        @test Gy.meshdim == 3
        # y-coords should contain the extrusion values
        @test minimum(Gy.nodes.coords[:, 2]) ≈ 0.0
        @test maximum(Gy.nodes.coords[:, 2]) ≈ 2.0

        # Extrude in z direction (default)
        Gz = extrude(G2d, 2.0, 1; dim=:z)
        @test Gz.meshdim == 3
        @test minimum(Gz.nodes.coords[:, 3]) ≈ 0.0
        @test maximum(Gz.nodes.coords[:, 3]) ≈ 2.0
    end

    @testset "Extrude - tag propagation" begin
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        G2d, _, _ = composite_voronoi_mesh_2d([0.1, 0.1], [1.0, 1.0];
            cell_constraints=wl, face_constraints=fl)

        G3d = extrude(G2d, 1.0, 2)

        # Cell tags should be propagated to all layers
        n_tagged_2d = count(G2d.cells.tag)
        n_tagged_3d = count(G3d.cells.tag)
        @test n_tagged_3d == n_tagged_2d * 2

        # Face tags should be propagated to vertical faces
        n_ftagged_2d = count(G2d.faces.tag)
        n_ftagged_3d = count(G3d.faces.tag)
        @test n_ftagged_3d == n_ftagged_2d * 2
    end

    @testset "Extrude - 3D geometry correctness" begin
        # Use composite mesh for a well-defined 2D mesh
        G2d, _, _ = composite_voronoi_mesh_2d([0.2, 0.2], [1.0, 1.0])
        compute_geometry!(G2d)

        G3d = extrude(G2d, 2.0, 1)
        compute_geometry!(G3d)

        # Total volume should be 2D area * height
        @test isapprox(sum(G3d.cells.volumes), sum(G2d.cells.volumes) * 2.0; atol=0.1)

        # All volumes should be positive
        @test all(G3d.cells.volumes .> 0)

        # Cell centroids z-coordinate should be at midpoint (1.0)
        @test all(isapprox.(G3d.cells.centroids[:, 3], 1.0; atol=1e-10))

        # All face areas should be positive
        @test all(G3d.faces.areas .> 0)
    end

    @testset "DistMesh 2D - basic" begin
        # Simple rectangular domain
        fd = p -> PolytopeMeshes._drectangle(p, 0.0, 1.0, 0.0, 1.0)
        fh = p -> ones(size(p, 1))
        bbox = [0.0 0.0; 1.0 1.0]
        pfix = zeros(0, 2)

        pts, tri = PolytopeMeshes.distmesh_2d(fd, fh, 0.2, bbox, pfix; max_iter=100)
        @test size(pts, 1) > 0
        @test size(pts, 2) == 2
        @test size(tri, 1) > 0
        @test size(tri, 2) == 3
    end

    @testset "DistMesh 2D - with fixed points" begin
        fd = p -> PolytopeMeshes._drectangle(p, 0.0, 1.0, 0.0, 1.0)
        fh = p -> ones(size(p, 1))
        bbox = [0.0 0.0; 1.0 1.0]
        pfix = [0.5 0.5; 0.25 0.25]

        pts, _ = PolytopeMeshes.distmesh_2d(fd, fh, 0.3, bbox, pfix; max_iter=100)
        # Fixed points should be the first rows
        @test isapprox(pts[1, :], [0.5, 0.5]; atol=1e-10)
        @test isapprox(pts[2, :], [0.25, 0.25]; atol=1e-10)
    end

    @testset "Distance functions" begin
        # drectangle
        p = [0.5 0.5; 1.5 0.5; 0.0 0.0]
        d = PolytopeMeshes._drectangle(p, 0.0, 1.0, 0.0, 1.0)
        @test d[1] < 0  # inside
        @test d[2] > 0  # outside

        # dpoly
        pv = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
        d2 = PolytopeMeshes._dpoly(p, pv)
        @test d2[1] < 0  # inside
        @test d2[2] > 0  # outside
    end

    @testset "Voronoi mesh 2D - basic" begin
        G, pts, F = voronoi_mesh_2d(0.2, [1.0, 1.0])

        @test G.cells.num > 0
        @test G.faces.num > 0
        @test G.nodes.num > 0

        compute_geometry!(G)
        @test isapprox(sum(G.cells.volumes), 1.0; atol=0.05)
        @test all(G.cells.volumes .> 0)
    end

    @testset "Voronoi mesh 2D - with cell constraints" begin
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = voronoi_mesh_2d(0.2, [1.0, 1.0];
            cell_constraints=wl)

        @test G.cells.num > 0
        @test any(G.cells.tag)  # Some cells should be tagged as wells
    end

    @testset "Voronoi mesh 2D - with face constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        G, pts, F = voronoi_mesh_2d(0.2, [1.0, 1.0];
            face_constraints=fl)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0
        @test any(G.faces.tag)  # Some faces should be tagged
    end

    @testset "Voronoi mesh 2D - with both constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = voronoi_mesh_2d(0.1, [1.0, 1.0];
            cell_constraints=wl, face_constraints=fl)

        @test G.cells.num > 0
        @test any(G.cells.tag)
        @test any(G.faces.tag)

        compute_geometry!(G)
        @test all(G.cells.volumes .> 0)
        @test isapprox(sum(G.cells.volumes), 1.0; atol=0.05)
    end

    @testset "Voronoi mesh 2D - polygon boundary" begin
        poly = Float64[0.0 0.0; 1.0 0.0; 0.5 1.0]
        G, _, _ = voronoi_mesh_2d(0.2, [1.0, 1.0]; poly_bdr=poly)

        @test G.cells.num > 0
        compute_geometry!(G)
        @test all(G.cells.volumes .> 0)
        @test isapprox(sum(G.cells.volumes), 0.5; atol=0.1)
    end

    @testset "Voronoi mesh 2D - refinement" begin
        wl = [Float64[0.5 0.2; 0.5 0.8]]
        G, _, _ = voronoi_mesh_2d(0.2, [1.0, 1.0];
            cell_constraints=wl, cc_refinement=true, cc_factor=0.5)

        @test G.cells.num > 0
        @test any(G.cells.tag)
    end

    @testset "Voronoi mesh 2D - suf_fc_cond=false" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        G, _, F = voronoi_mesh_2d(0.2, [1.0, 1.0];
            face_constraints=fl, suf_fc_cond=false)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0
    end

end
