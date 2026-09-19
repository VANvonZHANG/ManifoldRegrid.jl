using Aqua
using ManifoldRegrid
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(ManifoldRegrid)
end

@testset "ManifoldRegrid.jl" begin
    include("helpers.jl")
    include("test_schemes.jl")
end
