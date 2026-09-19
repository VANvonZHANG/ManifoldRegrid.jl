using ManifoldRegrid
using Test

@testset "scheme tags" begin
    @test Conservative() isa AbstractRemapScheme
    @test Bilinear() isa AbstractRemapScheme
end
