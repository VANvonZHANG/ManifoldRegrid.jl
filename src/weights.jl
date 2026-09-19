# Conservative remapping weights: overlap areas between source and destination
# cells, assembled into a sparse operator.

"""
    ConservativeWeights

First-order (piecewise-constant) conservative remapping operator between two
meshes covering the same sphere.

- `W`: sparse overlap-area matrix of size `(num_cells(dest), num_cells(src))`;
  `W[d, s]` is the area of `cell_src(s) ∩ cell_dest(d)`.
- `src`, `dest`: the meshes the weights were built from.
- `dst_areas`: `cell_volume(dest, d)` cache used for row normalization.

Construct via `conservative_weights`. The bare 4-field constructor performs no
validation, so it will happily build a `ConservativeWeights` that violates the
invariants below; `conservative_weights` is what checks them.

`W` always holds raw overlap areas so its invariants stay checkable: row sums
equal destination cell areas, column sums equal source cell areas. `remap`
applies the `1 / A_dest` normalization at application time.
"""
struct ConservativeWeights{W <: SparseMatrixCSC{Float64, Int},
    S <: AbstractManifoldMesh, D <: AbstractManifoldMesh}
    W::W
    src::S
    dest::D
    dst_areas::Vector{Float64}
end

# Angular radius of the smallest cap around the cell centroid containing the
# cell — the max over its vertices, which bounds the cell for convex cells
# lying within a hemisphere of the centroid (true for every shipped grid type;
# a violation surfaces as a conservation failure, not silently).
function _cell_circumradius(g::AbstractManifoldMesh, cell_id::Int)
    c = normalize(SVector{3, Float64}(cell_centroid(g, cell_id)))
    r = 0.0
    for nid in cell_nodes(g, cell_id)
        v = normalize(SVector{3, Float64}(node_coordinates(g, nid)))
        r = max(r, acos(clamp(dot(c, v), -1.0, 1.0)))
    end
    return r
end

function _check_conservation(w::ConservativeWeights)
    rtol = 1e-8
    rowsums = vec(sum(w.W; dims = 2))
    for d in eachindex(w.dst_areas)
        isapprox(rowsums[d], w.dst_areas[d]; rtol = rtol) || throw(ErrorException(
            "conservative weights are not conserved: destination cell $d " *
            "overlap sum $(rowsums[d]) != cell area $(w.dst_areas[d])"
        ))
    end
    colsums = vec(sum(w.W; dims = 1))
    for s in 1:num_cells(w.src)
        a = cell_volume(w.src, s)
        isapprox(colsums[s], a; rtol = rtol) || throw(ErrorException(
            "conservative weights are not conserved: source cell $s " *
            "overlap sum $(colsums[s]) != cell area $a"
        ))
    end
    return w
end

"""
    conservative_weights(src, dest) -> ConservativeWeights

Overlap-area weights for first-order conservative remapping from `src` to
`dest`. Both meshes must cover the same sphere with the same radius.

Destination cells are processed one at a time: candidate source cells come from
a k-d tree over source centroids with the circumcircle bound (never misses an
overlap), then each candidate is clipped exactly with
`spherical_polygon_intersection`. Construction verifies the conservation
invariants — row sums equal destination areas, column sums equal source areas —
to `rtol = 1e-8` and throws otherwise.
"""
function conservative_weights(src::AbstractManifoldMesh, dest::AbstractManifoldMesh)
    src.R == dest.R || throw(ArgumentError(
        "source and destination sphere radii differ: $(src.R) vs $(dest.R)"
    ))
    R = Float64(src.R)
    nc_src = num_cells(src)
    nc_dst = num_cells(dest)

    src_centroids = all_cell_centroids(src)
    dst_centroids = all_cell_centroids(dest)
    # Centroid scale depends on the grid type: the parametric grids return UNIT
    # vectors regardless of R, but `UnstructuredMesh` stores `R * P(cunit)` and
    # so returns R-scaled centroids. The `normalize` calls below are therefore
    # LOAD-BEARING, not cosmetic — they put the tree and the query point on the
    # unit sphere so the search radius is an angular chord. Deleting them breaks
    # remapping for `UnstructuredMesh` at R != 1.
    src_centroids = [normalize(SVector{3, Float64}(c)) for c in src_centroids]
    dst_centroids = [normalize(SVector{3, Float64}(c)) for c in dst_centroids]
    max_r_src = maximum(_cell_circumradius(src, c) for c in 1:nc_src)
    dst_r = [_cell_circumradius(dest, c) for c in 1:nc_dst]
    tree = KDTree(reduce(hcat, src_centroids))

    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for d in 1:nc_dst
        ring_d = cell_ring(dest, d)
        θ = min(dst_r[d] + max_r_src, π)
        radius = 2 * sin(θ / 2)   # unit-sphere chord; scale-free in R
        for s in inrange(tree, dst_centroids[d], radius)
            overlap = spherical_polygon_intersection(cell_ring(src, s), ring_d)
            area = spherical_polygon_area(overlap, R)
            if area > 0
                push!(rows, d)
                push!(cols, s)
                push!(vals, area)
            end
        end
    end

    W = sparse(rows, cols, vals, nc_dst, nc_src)
    w = ConservativeWeights(W, src, dest, [cell_volume(dest, d) for d in 1:nc_dst])
    return _check_conservation(w)
end
