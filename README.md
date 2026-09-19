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

## Conservative remapping from a model grid

Unstructured model output (ICON, MPAS, …) reads through `ManifoldFields` and
remaps without any per-grid code:

```julia
using ManifoldFields, ManifoldMeshes, ManifoldRegrid

fs = load_ugrid("oQU480.ugrid.nc")            # FieldSet on an UnstructuredMesh
dst = LatLonGrid(lat_edges = collect(-90.0:1.0:90.0),
                 lon_edges = collect(0.0:1.0:360.0))

w = conservative_weights(mesh(fs), dst)       # build once
out = remap(w, fs[:bottomDepth])              # DiscreteField{CellLoc} on dst
```

Global integrals of `fs[:bottomDepth]` and `out` agree to machine precision, as
long as the two meshes cover the same region — `conservative_weights` verifies
that at construction. The caveat bites for this particular file: the tutorial's
`oQU480.ugrid.nc` is an ocean-only mesh (`Σ cell_volume ≈ 0.70 × 4π`), so a
global destination grid is rejected and `dst` must be chosen over the source
footprint. The example file is the one shipped with the `ManifoldFields.jl`
tutorial data (`ManifoldFields.jl/examples/data/oQU480.ugrid.nc`), so this
recipe is not part of the test suite — the suite uses a self-contained
octahedral `UnstructuredMesh` instead.

**Package status:** private/development package; not registered. See `ROADMAP.md`.
