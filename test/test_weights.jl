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

@testset "radius mismatch" begin
    g = coarse_grid()
    g2 = LatLonGrid(
        lat_edges = collect(-90.0:30.0:90.0),
        lon_edges = collect(0.0:45.0:360.0),
        R = 2.0
    )
    @test_throws ArgumentError conservative_weights(g, g2)
end
