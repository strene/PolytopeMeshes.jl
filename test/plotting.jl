using CairoMakie

@testset "Plotting (Makie extension)" begin
    @test PolytopeMeshes.check_plotting_availability()
    @test !PolytopeMeshes.plotting_is_interactive()
    fl = [[[0.2, 0.2], [0.8, 0.8]]]
    wl = [[[0.2, 0.8], [0.8, 0.2]]]
    m2 = composite_pebi_grid_2d([0.1, 0.1], [1, 1]; cell_constraints = wl, face_constraints = fl)
    m3 = composite_pebi_grid_3d((4, 4, 4), (1, 1, 1); cell_constraints = [[[0.5, 0.5, 0.1], [0.5, 0.5, 0.9]]])
    mktempdir() do dir
        for (name, m) in (("m2", m2), ("m3", m3))
            fig, ax, p = plot_mesh(m)
            @test fig isa Figure
            save(joinpath(dir, "$name.png"), fig)
            @test isfile(joinpath(dir, "$name.png"))
            fig, ax, p = plot_cell_data(m, [s[1] for s in m.cell_sites]; cells = 1:3)
            @test fig isa Figure
            fig = Figure()
            ax = m isa PolytopalMesh{2} ? Axis(fig[1, 1]) : Axis3(fig[1, 1])
            @test !isnothing(plot_mesh!(ax, m; edges = false, constraints = false))
            @test isnothing(plot_faces!(ax, m; faces = Int[]))
            fig = @test_logs (:warn,) match_mode = :any plot_interactive(m; data = Dict("x" => [s[1] for s in m.cell_sites]))
            @test fig isa Figure
            # Drive the widgets: select a subset, an empty set and another field
            for s in [c for c in fig.content if c isa IntervalSlider]
                set_close_to!(s, 0.0, 0.5)
            end
            for t in [c for c in fig.content if c isa Toggle]
                t.active[] = !t.active[]
            end
            menu = only([c for c in fig.content if c isa Menu])
            menu.i_selected[] = length(menu.options[])
            save(joinpath(dir, "$(name)_interactive.png"), fig)
            @test isfile(joinpath(dir, "$(name)_interactive.png"))
        end
        @test_throws ArgumentError plot_interactive(m2; data = Dict("bad" => [1.0]))
    end
end
