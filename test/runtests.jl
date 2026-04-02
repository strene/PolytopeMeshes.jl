using Test
using PolytopeMeshes

@testset "PolytopeMeshes.jl" begin

    @testset "Grid data structure" begin
        # Create a simple 2-cell grid (two triangles)
        coords = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.5]
        cell_faces = [1, 2, 3, 4, 5, 3]
        cell_facePos = [1, 4, 7]
        face_nodes = [1, 2, 2, 5, 5, 1, 2, 3, 3, 5]
        face_nodePos = [1, 3, 5, 7, 9, 11]
        face_neighbors = [1 0; 1 2; 1 0; 2 0; 2 0]

        G = UnstructuredGrid(coords, cell_faces, cell_facePos, face_nodes,
            face_nodePos, face_neighbors; griddim=2)

        @test G.cells.num == 2
        @test G.faces.num == 5
        @test G.nodes.num == 5
        @test G.griddim == 2
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

    @testset "Clipped PEBI 2D - simple" begin
        # Regular grid of points
        pts = zeros(9, 2)
        idx = 1
        for y in [0.25, 0.5, 0.75]
            for x in [0.25, 0.5, 0.75]
                pts[idx, :] = [x, y]
                idx += 1
            end
        end
        bnd = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]

        G = clipped_pebi_2d(pts, bnd)

        @test G.cells.num == 9
        @test G.faces.num > 0
        @test G.nodes.num > 0
        @test G.griddim == 2

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

    @testset "Composite PEBI grid 2D - basic" begin
        G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0])

        @test G.cells.num > 0
        @test G.faces.num > 0
        @test G.nodes.num > 0

        compute_geometry!(G)
        @test isapprox(sum(G.cells.volumes), 1.0; atol=0.1)
    end

    @testset "Composite PEBI grid 2D - with cell constraints" begin
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
            cell_constraints=wl)

        @test G.cells.num > 0
        @test any(G.cells.tag)  # Some cells should be tagged as wells
    end

    @testset "Composite PEBI grid 2D - with face constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
            face_constraints=fl)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0
    end

    @testset "Composite PEBI grid 2D - with both constraints" begin
        fl = [Float64[0.2 0.2; 0.8 0.8]]
        wl = [Float64[0.2 0.8; 0.8 0.2]]
        G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
            cell_constraints=wl,
            face_constraints=fl)

        @test G.cells.num > 0
        @test size(F.f_pts, 1) > 0

        compute_geometry!(G)
        @test all(G.cells.volumes .> 0)
    end

end
