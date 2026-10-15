# Waterflooding on a 3D PEBI grid with JutulDarcy. The grid conforms to two
# vertical wells (cell constraints) and a vertical fault between them (face
# constraint). The fault is made sealing by setting the transmissibility of
# the fault faces to zero, which forces the water around the fault.
#
# Requires JutulDarcy (and GLMakie or CairoMakie for the plots), which are not
# dependencies of PolytopeMeshes:
#   using Pkg; Pkg.add(["JutulDarcy", "GLMakie"])
using PolytopeMeshes, Jutul, JutulDarcy
using GLMakie
Darcy, bar, day = si_units(:darcy, :bar, :day)

## Grid
L, H = 1000.0, 100.0
wells = [
    [[175.0, 175.0, 5.0], [175.0, 175.0, H - 5.0]], # Injector
    [[825.0, 825.0, 5.0], [825.0, 825.0, H - 5.0]], # Producer
]
fault = surface_triangulation([[250.0, 750.0, 0.0], [750.0, 250.0, 0.0], [750.0, 250.0, H], [250.0, 750.0, H]], 25.0)
G = composite_pebi_grid_3d((20, 20, 4), (L, L, H);
    cell_constraints = wells, cc_rho = H / 4,
    face_constraints = [fault], fc_rho = 25.0)

# Convert to a Jutul mesh. The cell and face constraints become mesh entity
# tags, which give the well cells and the fault faces.
g = UnstructuredMesh(G)
injector_cells = Jutul.get_mesh_entity_tag(g, Cells(), :cell_constraints, :constraint_1)
producer_cells = Jutul.get_mesh_entity_tag(g, Cells(), :cell_constraints, :constraint_2)
fault_faces = Jutul.get_mesh_entity_tag(g, Faces(), :face_constraints, :constraint_1)

## Reservoir model
domain = reservoir_domain(g, permeability = 0.5Darcy, porosity = 0.2)
Inj = setup_well(domain, injector_cells, name = :Injector)
Prod = setup_well(domain, producer_cells, name = :Producer)
rho = [1000.0, 800.0]
sys = ImmiscibleSystem((AqueousPhase(), LiquidPhase()), reference_densities = rho)
model, parameters = setup_reservoir_model(domain, sys, wells = [Inj, Prod], extra_out = true)
replace_variables!(model,
    RelativePermeabilities = BrooksCoreyRelativePermeabilities(sys, [2.0, 2.0]),
    PhaseMassDensities = ConstantCompressibilityDensities(
        p_ref = 100bar, density_ref = rho, compressibility = [1e-6, 1e-5] / bar)
)

dt = fill(30.0day, 40)
pv = pore_volume(model, parameters)
controls = Dict(
    :Injector => InjectorControl(TotalRateTarget(0.8 * sum(pv) / sum(dt)), [1.0, 0.0], density = rho[1]),
    :Producer => ProducerControl(BottomHolePressureTarget(80bar)),
)
forces = setup_reservoir_forces(model, control = controls)
state0 = setup_reservoir_state(model, Pressure = 100bar, Saturations = [0.0, 1.0])

## Simulate with an open and a sealing fault
open_fault = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces)
sealed = deepcopy(parameters)
sealed[:Reservoir][:Transmissibilities][fault_faces] .= 0.0
sealed_fault = simulate_reservoir(state0, model, dt, parameters = sealed, forces = forces)

## Results
fig = Figure(size = (900, 400))
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Water cut")
for (label, res) in (("Open fault", open_fault), ("Sealing fault", sealed_fault))
    lines!(ax, res.time ./ day, res.wells[:Producer][:wcut], label = label)
end
axislegend(ax, position = :lt)
display(fig)

# Final water saturation with the sealing fault, in the bottom half of the
# grid, together with the fault and the well cells
sw = sealed_fault.states[end][:Saturations][1, :]
fig, ax, plt = plot_polymesh_data(G, sw; cells = [s[3] < H / 2 for s in G.cell_sites], colormap = :Blues,
    axis = (aspect = (1, 1, 0.3),))
display(fig)
plot_polymesh_interactive(G; data = Dict("Sw (sealing fault)" => sw,
    "Sw (open fault)" => open_fault.states[end][:Saturations][1, :]))
