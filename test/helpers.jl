using DimensionalData
using ManifoldFields
using ManifoldMeshes

# Octahedron: 8 spherical triangles, each of area π/2 — the smallest
# UnstructuredMesh that tiles the sphere exactly.
function octahedron()
    UnstructuredMesh(
        [0.0, 90.0, 180.0, 270.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 90.0, -90.0],
        [1 2 5; 2 3 5; 3 4 5; 4 1 5; 2 1 6; 3 2 6; 4 3 6; 1 4 6];
        start_index = 1
    )
end

function coarse_grid()
    LatLonGrid(
        lat_edges = collect(-90.0:30.0:90.0),
        lon_edges = collect(0.0:45.0:360.0)
    )
end

function fine_grid()
    LatLonGrid(
        lat_edges = collect(-90.0:15.0:90.0),
        lon_edges = collect(0.0:15.0:360.0)
    )
end

function small_other_grid()
    LatLonGrid(
        lat_edges = collect(-90.0:45.0:90.0),
        lon_edges = collect(0.0:90.0:360.0)
    )
end
