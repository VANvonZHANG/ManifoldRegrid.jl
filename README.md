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

**Package status:** private/development package; not registered. See `ROADMAP.md`.
