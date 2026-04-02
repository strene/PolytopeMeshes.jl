# PolytopeMeshes.jl

A Julia package for generating conforming Voronoi meshes (PEBI grids) that adapt
to geological features such as faults, fractures, and wells.

Based on the UPR (Unstructured PEBI-grids for Reservoir) module from the
[MATLAB Reservoir Simulation Toolbox (MRST)](https://github.com/SINTEF-AppliedCompSci/MRST/tree/main/modules/upr).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/strene/PolytopeMeshes.jl")
```

## Quick Start

```julia
using PolytopeMeshes

# Create a basic Voronoi grid
G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0])

# Create a grid with a well (cell constraint) and a fault (face constraint)
well = [Float64[0.2 0.8; 0.8 0.2]]
fault = [Float64[0.2 0.2; 0.8 0.8]]
G, pts, F = composite_pebi_grid_2d([0.1, 0.1], [1.0, 1.0];
    cell_constraints=well,
    face_constraints=fault)

# Compute geometric quantities
compute_geometry!(G)

# Access grid data (MRST-compatible format)
G.cells.num        # Number of cells
G.cells.volumes    # Cell volumes
G.cells.centroids  # Cell centroids
G.cells.tag        # Tagged well cells
G.faces.num        # Number of faces
G.faces.neighbors  # Face-to-cell connectivity
G.nodes.coords     # Node coordinates
```

## Grid Data Structure

The grid follows the MRST unstructured grid format:

| Field | Description |
|-------|-------------|
| `G.cells.num` | Number of cells |
| `G.cells.faces` | Cell-to-face mapping |
| `G.cells.facePos` | Position pointers for cell-face map |
| `G.cells.volumes` | Cell volumes (after `compute_geometry!`) |
| `G.cells.centroids` | Cell centroids (after `compute_geometry!`) |
| `G.cells.tag` | Boolean tags for constraint cells |
| `G.faces.num` | Number of faces |
| `G.faces.nodes` | Face-to-node mapping |
| `G.faces.nodePos` | Position pointers for face-node map |
| `G.faces.neighbors` | Face neighbor cells (0 = boundary) |
| `G.faces.areas` | Face areas (after `compute_geometry!`) |
| `G.faces.normals` | Face normal vectors (after `compute_geometry!`) |
| `G.faces.tag` | Boolean tags for constraint faces |
| `G.nodes.coords` | Node coordinates [n × dim] |
| `G.nodes.num` | Number of nodes |

## Main Functions

### Grid Generation

- **`composite_pebi_grid_2d(celldim, pdims; kwargs...)`** — Main entry point.
  Creates a composite PEBI grid with a Cartesian background refined around constraints.

- **`clipped_pebi_2d(pts, bnd)`** — Creates a Voronoi grid from sites clipped to a polygon boundary.

### Constraint Site Generation

- **`line_sites_2d(constraints, grid_size; kwargs...)`** — Places Voronoi sites
  along cell constraints (wells). Sites are placed directly on the constraint
  path so that cell centroids trace the well.

- **`surface_sites_2d(constraints, grid_size; kwargs...)`** — Places Voronoi sites
  on both sides of face constraints (faults/fractures). Uses the circle-based
  bisection method to ensure fault faces align with the constraint.

- **`surface_suf_cond_2d(pts, F)`** — Enforces the sufficient surface condition
  by removing reservoir sites that lie inside fault circles.

### Utilities

- **`remove_conflict_points(P1, P2, dist)`** — Removes sites from P1 too close to P2.
- **`compute_geometry!(G)`** — Computes face areas/normals and cell volumes/centroids.

## References

- Berge, R.L., Klemetsdal, Ø.S. & Lie, K.-A. (2019). "Unstructured Voronoi grids
  conforming to lower dimensional objects." *Computational Geosciences*, 23, 169–188.
  https://doi.org/10.1007/s10596-018-9790-0

## License

See [LICENSE](LICENSE) for details.
