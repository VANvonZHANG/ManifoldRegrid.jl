using LinearAlgebra
using ManifoldRegrid
using ManifoldMeshes
using SparseArrays
using Test

@testset "conservative_weights invariants" begin
    src = coarse_grid()
    dst = fine_grid()
    w = conservative_weights(src, dst)

    @test size(w.W) == (num_cells(dst), num_cells(src))
    @test w.src === src
    @test w.dest === dst
    @test all(≥(0), nonzeros(w.W))

    rowsums = vec(sum(w.W; dims = 2))
    colsums = vec(sum(w.W; dims = 1))
    @test rowsums ≈ [cell_volume(dst, d) for d in 1:num_cells(dst)] rtol = 1e-8
    @test colsums ≈ [cell_volume(src, s) for s in 1:num_cells(src)] rtol = 1e-8
    @test sum(rowsums) ≈ 4π rtol = 1e-8

    # the checker must actually throw: the spec calls loud failure the product's
    # core promise, and nothing else exercises the throwing branch
    @test_throws ErrorException ManifoldRegrid._check_conservation(
        ConservativeWeights(sparse([1], [1], [1.0], size(w.W)...), w.src, w.dest,
        w.dst_areas)
    )
    # Clipping noise must not leak in as microscopic weights: every stored entry
    # is a real geometric sliver, so the bar is strictly positive rather than a
    # magnitude threshold. Sliver scale is a property of the fixture, not of the
    # code: the smallest one here is ~1.95e-3, and a finer destination — LatLon
    # -> CubedSphereGrid(n = 8) — drops it to ~3.06e-5. A threshold like 1e-6
    # would be a fixture-scale accident waiting to trip on a legitimate sliver.
    @test minimum(nonzeros(w.W)) > 0.0
end

# The suite's only non-LatLon fixtures are the octahedron (unstructured); the
# parametric grids other than LatLon are otherwise never a remapping partner.
# Pin the "works for every grid type" claim for them, in both roles.
@testset "CubedSphere and ReducedGaussian as source and destination" begin
    for g in (CubedSphereGrid(n = 4), ReducedGaussianGrid(nlat = 6))
        for (src, dst) in ((g, coarse_grid()), (coarse_grid(), g), (g, g))
            w = conservative_weights(src, dst)
            rowsums = vec(sum(w.W; dims = 2))
            colsums = vec(sum(w.W; dims = 1))
            @test rowsums ≈ [cell_volume(dst, d) for d in 1:num_cells(dst)] rtol = 1e-8
            @test colsums ≈ [cell_volume(src, s) for s in 1:num_cells(src)] rtol = 1e-8
        end
    end
end

# Clipping — unlike the area formulas — depends on the ring orientation, so the
# counter-clockwise guarantee gets checked directly: the cell centroid must be
# inside every edge of its own ring.
@testset "cell rings are counter-clockwise ($(nameof(typeof(g))))" for g in (
    CubedSphereGrid(n = 4), ReducedGaussianGrid(nlat = 6))
    for c in 1:num_cells(g)
        ring = cell_ring(g, c)
        centroid = normalize(cell_centroid(g, c))
        inside = true
        for i in eachindex(ring)
            a = ring[i]
            b = ring[mod1(i + 1, length(ring))]
            inside &= side_of_geodesic(centroid, a, b) >= 0
        end
        @test inside
    end
end

@testset "same-mesh weights are diagonal" begin
    g = coarse_grid()
    w = conservative_weights(g, g)
    @test nnz(w.W) == num_cells(g)
    @test diag(w.W) ≈ [cell_volume(g, c) for c in 1:num_cells(g)] rtol = 1e-12
end

@testset "unstructured source (octahedron)" begin
    oct = octahedron()
    dst = coarse_grid()
    w = conservative_weights(oct, dst)

    @test vec(sum(w.W; dims = 1)) ≈ fill(π / 2, num_cells(oct)) rtol = 1e-8
    @test vec(sum(w.W; dims = 2)) ≈ [cell_volume(dst, d) for d in 1:num_cells(dst)] rtol = 1e-8
end

@testset "unstructured destination (round trip)" begin
    oct = octahedron()
    w = conservative_weights(coarse_grid(), oct)
    @test size(w.W) == (num_cells(oct), num_cells(coarse_grid()))
    @test vec(sum(w.W; dims = 2)) ≈ fill(π / 2, num_cells(oct)) rtol = 1e-8
end

@testset "non-unit radius (R != 1)" begin
    # R != 1 is the discriminating choice: the candidate-search radius used to
    # be `2R * sin(θ / 2)`, which is the unit-sphere chord only when R == 1.
    R = 0.5
    src = LatLonGrid(
        lat_edges = collect(-90.0:30.0:90.0),
        lon_edges = collect(0.0:45.0:360.0),
        R = R
    )
    dst = LatLonGrid(
        lat_edges = collect(-90.0:15.0:90.0),
        lon_edges = collect(0.0:15.0:360.0),
        R = R
    )
    w = conservative_weights(src, dst)

    rowsums = vec(sum(w.W; dims = 2))
    colsums = vec(sum(w.W; dims = 1))
    @test rowsums ≈ [cell_volume(dst, d) for d in 1:num_cells(dst)] rtol = 1e-8
    @test colsums ≈ [cell_volume(src, s) for s in 1:num_cells(src)] rtol = 1e-8
    @test sum(rowsums) ≈ 4π * R^2 rtol = 1e-8
end

@testset "radius mismatch" begin
    g = coarse_grid()
    g2 = LatLonGrid(
        lat_edges = collect(-90.0:30.0:90.0),
        lon_edges = collect(0.0:45.0:360.0),
        R = 2.0
    )
    @test_throws ArgumentError conservative_weights(g, g2)
end
