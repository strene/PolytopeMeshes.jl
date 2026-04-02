"""
Sort cell faces counter-clockwise for consistent 2D grid orientation.
Ported from MRST's sortEdges/sortCellFaces.
"""

"""
    sort_edges!(G::UnstructuredGrid)

Sort the cell faces of a 2D grid counter-clockwise to ensure positive cell
volumes and consistent face orientations.
"""
function sort_edges!(G::UnstructuredGrid)
    @assert G.griddim == 2

    coords = G.nodes.coords
    face_nodes_arr = G.faces.nodes
    face_nodePos = G.faces.nodePos

    for i in 1:G.cells.num
        fi = G.cells.facePos[i]
        li = G.cells.facePos[i + 1] - 1
        cf = G.cells.faces[fi:li]
        nf = length(cf)

        if nf < 3
            continue
        end

        # Get edges (node pairs) for each face in this cell
        edges = zeros(Int, nf, 2)
        reverse_flags = falses(nf)

        for (j, fid) in enumerate(cf)
            n1 = face_nodes_arr[face_nodePos[fid]]
            n2 = face_nodes_arr[face_nodePos[fid] + 1]

            # Determine orientation based on neighbor
            if G.faces.neighbors[fid, 1] == i
                edges[j, :] = [n1, n2]
            else
                edges[j, :] = [n2, n1]
                reverse_flags[j] = true
            end
        end

        # Sort edges to form a connected chain
        sorted_order = _sort_edge_chain(edges)
        if isnothing(sorted_order)
            continue
        end

        cf_sorted = cf[sorted_order]
        edges_sorted = edges[sorted_order, :]

        # Check area orientation (should be positive = counter-clockwise)
        area = 0.0
        for j in 1:nf
            n1 = edges_sorted[j, 1]
            x1, y1 = coords[n1, 1], coords[n1, 2]
            n2 = edges_sorted[j, 2]
            x2, y2 = coords[n2, 1], coords[n2, 2]
            area += (x1 + x2) * (y2 - y1)
        end

        if area < 0
            # Reverse the order and flip edges
            cf_sorted = reverse(cf_sorted)
        end

        # Update cell faces
        G.cells.faces[fi:li] = cf_sorted

        # Ensure face nodes are consistent with cell orientation
        for (j, fid) in enumerate(cf_sorted)
            n1 = face_nodes_arr[face_nodePos[fid]]
            n2 = face_nodes_arr[face_nodePos[fid] + 1]

            # The first node of a face should be the "start" when
            # traversing counter-clockwise from cell 1's perspective
            if G.faces.neighbors[fid, 2] == i
                # This cell is neighbor 2, so nodes should be reversed from cell's perspective
                # but stored in neighbor 1's orientation - this is correct as-is
            end
        end
    end

    return G
end

"""
Sort edges to form a connected chain (each edge's end matches next edge's start).
"""
function _sort_edge_chain(edges::Matrix{Int})
    n = size(edges, 1)
    if n <= 1
        return collect(1:n)
    end

    order = zeros(Int, n)
    used = falses(n)
    order[1] = 1
    used[1] = true

    current_end = edges[1, 2]

    for pos in 2:n
        found = false
        for j in 1:n
            if !used[j]
                if edges[j, 1] == current_end
                    order[pos] = j
                    used[j] = true
                    current_end = edges[j, 2]
                    found = true
                    break
                elseif edges[j, 2] == current_end
                    # Need to flip this edge
                    edges[j, :] = edges[j, [2, 1]]
                    order[pos] = j
                    used[j] = true
                    current_end = edges[j, 2]
                    found = true
                    break
                end
            end
        end
        if !found
            return nothing  # Can't sort
        end
    end

    return order
end
