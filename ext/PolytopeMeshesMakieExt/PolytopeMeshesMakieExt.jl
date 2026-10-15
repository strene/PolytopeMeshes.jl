module PolytopeMeshesMakieExt

using PolytopeMeshes
using PolytopeMeshes: PolytopalMesh, plot_triangulation, number_of_cells, number_of_faces, face_nodes, cell_faces
using Makie
using StaticArrays

PolytopeMeshes.check_plotting_availability_impl() = true

function PolytopeMeshes.current_backend_name()
    backend = Makie.current_backend()
    return isnothing(backend) || backend isa Missing ? :none : nameof(backend)
end

include("mesh_plots.jl")
include("interactive.jl")

end # module
