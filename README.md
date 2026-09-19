# ManifoldRegrid.jl

Conservative and bilinear remapping between `ManifoldMeshes.jl` meshes.

```julia
using ManifoldFields
using ManifoldMeshes
using ManifoldRegrid
using DimensionalData

src = LatLonGrid(lat_edges = collect(-90.0:30.0:90.0), lon_edges = collect(0.0:45.0:360.0))
dst = LatLonGrid(lat_edges = collect(-90.0:15.0:90.0), lon_edges = collect(0.0:15.0:360.0))

values = randn(num_cells(src))
f = DiscreteField(CellLoc, src, values, (Dim{:cell}(1:num_cells(src)),); name = :emission)

coarse_to_fine = remap(Conservative(), f, dst)   # DiscreteField{CellLoc} on dst
```

## API

```julia
w = conservative_weights(src, dst)   # build once, reuse across fields and time steps
remap(w, f)                          # apply to a DiscreteField{CellLoc}
remap(Conservative(), f, dst)        # one-shot: builds weights then applies
remap(Bilinear(), f, dst)            # nodal bilinear (== ManifoldFields.interpolate(f, dst))
```

`Conservative` only accepts `DiscreteField{CellLoc}` and `Bilinear` only
`DiscreteField{NodeLoc}`; any other combination is a `MethodError`.

Unstructured sources work unchanged — e.g. an ICON/MPAS-style mesh loaded
through `ManifoldFields.load_ugrid` reconstructs as a
`ManifoldMeshes.UnstructuredMesh`, which is a first-class source mesh here.

**Package status:** private/development package; not registered. See `ROADMAP.md`.
