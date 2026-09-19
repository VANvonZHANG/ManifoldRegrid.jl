module ManifoldRegrid

using DimensionalData
using LinearAlgebra
using ManifoldFields
using ManifoldMeshes
using NearestNeighbors
using SparseArrays
using StaticArrays

export AbstractRemapScheme, Conservative, Bilinear
export ConservativeWeights, conservative_weights
export remap

"""
    AbstractRemapScheme

Supertype of remapping schemes. A scheme is a dispatch tag: `remap(scheme, f,
dest_mesh)` is only defined for the field locations the scheme supports, so an
unsupported combination is a `MethodError` at the call site.
"""
abstract type AbstractRemapScheme end

"""
    Conservative()

First-order (piecewise-constant) conservative remapping. Applies to
cell-centred fields; destination values are area-weighted means of the source
cells overlapping them, and the global integral is preserved.
"""
struct Conservative <: AbstractRemapScheme end

"""
    Bilinear()

Nodal bilinear remapping: destination node values are bilinear interpolations
of the source nodal field.
"""
struct Bilinear <: AbstractRemapScheme end

include("weights.jl")
include("remap.jl")

end # module
