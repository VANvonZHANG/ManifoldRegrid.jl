# Remapping front-end. Schemes dispatch on the field's location tag, so an
# unsupported combination is a MethodError at the call site rather than a
# runtime surprise.

"""
    remap(scheme, f, dest_mesh) -> DiscreteField

Remap `f` onto `dest_mesh` with `scheme`.

- `remap(Bilinear(), f::DiscreteField{NodeLoc}, dest)` — nodal bilinear
  interpolation, equivalent to `ManifoldFields.interpolate(f, dest)`.
- `remap(Conservative(), f::DiscreteField{CellLoc}, dest)` — first-order
  conservative remapping from `conservative_weights(mesh(f), dest)`.

The scheme and the field location must match: `Conservative` is only defined
for cell-centred fields and `Bilinear` only for node-centred ones, so any other
combination raises a `MethodError`.
"""
function remap end

function remap(::Bilinear, f::DiscreteField{NodeLoc}, dest::AbstractManifoldMesh)
    return interpolate(f, dest)
end

function remap(::Conservative, f::DiscreteField{CellLoc}, dest::AbstractManifoldMesh)
    return remap(conservative_weights(mesh(f), dest), f)
end

"""
    remap(w::ConservativeWeights, f::DiscreteField{CellLoc}) -> DiscreteField{CellLoc}

Apply precomputed conservative weights: each destination value is the
area-weighted mean of the source values overlapping it,

    dest[d] = Σₛ W[d, s] * src[s] / cell_volume(dest, d).

Dividing by the exact destination cell area rather than by the row sum of `W`
keeps the global integral conserved to machine precision even when the overlap
areas carry rounding error: geometric error degrades local accuracy only.

The field must be defined on the *same mesh object* the weights were built
from (`mesh(f) === w.src`); an independently reconstructed but geometrically
identical mesh is rejected — rebuild the weights alongside it.

The result is always cell-major: the location axis comes first and trailing
axes keep their relative order, so a `(time, cell)` field comes back
`(cell, time)`.
"""
function remap(w::ConservativeWeights, f::DiscreteField{CellLoc})
    mesh(f) === w.src || throw(ArgumentError(
        "field is not defined on the weights' source mesh"
    ))
    values = data(f)
    loc_axis = location_axis(f)
    trailing_dims = Tuple(d for (i, d) in enumerate(dims(f)) if i != loc_axis)
    other_axes = Tuple(i for i in 1:ndims(values) if i != loc_axis)

    moved = permutedims(values, (loc_axis, other_axes...))
    flat = reshape(moved, num_cells(w.src), :)
    out = reshape((w.W * flat) ./ w.dst_areas,
        (num_cells(w.dest), map(length, trailing_dims)...))

    cell_dim = location_dimname(CellLoc)
    return DiscreteField(CellLoc, w.dest, out,
        (cell_dim(1:num_cells(w.dest)), trailing_dims...);
        name = DimensionalData.name(f))
end
