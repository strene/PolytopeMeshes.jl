module PolytopeMeshesJutulExt

using PolytopeMeshes
import Jutul

"""
    jutul_mesh(G::UnstructuredMesh)

Convert a PolytopeMeshes `UnstructuredMesh` to a Jutul `UnstructuredMesh`.

Only 3D meshes are supported by Jutul. If a 2D mesh is provided, it is
automatically extruded to 3D with a single layer of length 1 in the
z-direction (equivalent to `extrude(G, 1.0, 1)`).
"""
function PolytopeMeshes.jutul_mesh(G::UnstructuredMesh)
    if G.meshdim == 2
        G = extrude(G, 1.0, 1)
    end

    # Jutul's outer constructor expects (2 × nf) for neighbors
    # and (dim × nn) for node coordinates
    N = collect(G.faces.neighbors')
    coords = collect(G.nodes.coords')

    return Jutul.UnstructuredMesh(
        G.cells.faces,
        G.cells.facePos,
        G.faces.nodes,
        G.faces.nodePos,
        coords,
        N
    )
end

end # module PolytopeMeshesJutulExt
