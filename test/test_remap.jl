using DimensionalData
using ManifoldFields
using ManifoldMeshes
using ManifoldRegrid
using Test

@testset "conservative remap reproduces constants" begin
    src = coarse_grid()
    dst = fine_grid()
    f = DiscreteField(
        CellLoc, src, fill(3.0, num_cells(src)), (Dim{:cell}(1:num_cells(src)),);
        name = :emission
    )
    out = remap(Conservative(), f, dst)

    @test out isa DiscreteField{CellLoc}
    @test mesh(out) === dst
    @test DimensionalData.name(out) == :emission
    @test data(out) ≈ fill(3.0, num_cells(dst))
end

@testset "global integral is conserved" begin
    src = coarse_grid()
    dst = fine_grid()
    # The bare pattern integrates to exactly zero over the sphere, which leaves
    # both totals as roundoff noise (~1e-17) and makes an rtol comparison
    # meaningless; the +1 gives the check a nonnegligible reference (4π).
    latlons = [ManifoldMeshes._cartesian_to_latlon(cell_centroid(src, c))
               for c in 1:num_cells(src)]
    values = [1.0 + sin(deg2rad(lat)) * cos(deg2rad(lon)) for (lat, lon) in latlons]
    f = DiscreteField(CellLoc, src, values, (Dim{:cell}(1:num_cells(src)),); name = :x)

    out = remap(Conservative(), f, dst)
    src_total = sum(values .* [cell_volume(src, s) for s in 1:num_cells(src)])
    dst_total = sum(data(out) .* [cell_volume(dst, d) for d in 1:num_cells(dst)])
    @test dst_total ≈ src_total rtol = 1e-10
end

@testset "trailing dimensions" begin
    src = coarse_grid()
    dst = fine_grid()
    n = num_cells(src)
    values = reshape(collect(1.0:(n * 3)), n, 3)
    ft = DiscreteField(CellLoc, src, values,
        (Dim{:cell}(1:n), Dim{:time}(1:3)); name = :x)
    outt = remap(Conservative(), ft, dst)

    @test size(data(outt)) == (num_cells(dst), 3)
    for t in 1:3
        col = DiscreteField(CellLoc, src, values[:, t], (Dim{:cell}(1:n),); name = :x)
        @test data(outt)[:, t] ≈ data(remap(Conservative(), col, dst))
    end
end

@testset "reusing precomputed weights" begin
    src = coarse_grid()
    dst = fine_grid()
    w = conservative_weights(src, dst)
    f = DiscreteField(CellLoc, src, fill(7.0, num_cells(src)),
        (Dim{:cell}(1:num_cells(src)),); name = :x)

    @test data(remap(w, f)) ≈ data(remap(Conservative(), f, dst))

    other = small_other_grid()
    g = DiscreteField(CellLoc, other, fill(1.0, num_cells(other)),
        (Dim{:cell}(1:num_cells(other)),); name = :x)
    @test_throws ArgumentError remap(w, g)
end

@testset "bilinear remap matches interpolate" begin
    src = coarse_grid()
    dst = fine_grid()
    f = DiscreteField(NodeLoc, src, collect(1.0:num_nodes(src)),
        (Dim{:node}(1:num_nodes(src)),); name = :x)

    out = remap(Bilinear(), f, dst)
    @test out isa DiscreteField{NodeLoc}
    @test data(out) ≈ data(interpolate(f, dst))
end

@testset "scheme / location mismatches are MethodErrors" begin
    src = coarse_grid()
    dst = fine_grid()
    node_field = DiscreteField(NodeLoc, src, collect(1.0:num_nodes(src)),
        (Dim{:node}(1:num_nodes(src)),); name = :x)
    cell_field = DiscreteField(CellLoc, src, collect(1.0:num_cells(src)),
        (Dim{:cell}(1:num_cells(src)),); name = :x)

    @test_throws MethodError remap(Conservative(), node_field, dst)
    @test_throws MethodError remap(Bilinear(), cell_field, dst)
end
