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

## Conservative remapping from an unstructured source

Unstructured model output (ICON, MPAS, …) reads through `ManifoldFields` and
remaps without any per-grid code. The runnable form of that — a *global*
unstructured source — is the octahedral mesh the test suite builds:

```julia
using DimensionalData
using ManifoldFields, ManifoldMeshes, ManifoldRegrid

# 8 spherical triangles, Σ cell_volume = 4π: a global UnstructuredMesh
oct = UnstructuredMesh(
    [0.0, 90.0, 180.0, 270.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 90.0, -90.0],
    [1 2 5; 2 3 5; 3 4 5; 4 1 5; 2 1 6; 3 2 6; 4 3 6; 1 4 6];
    start_index = 1
)
f = DiscreteField(CellLoc, oct, fill(1.0, num_cells(oct)),
    (Dim{:cell}(1:num_cells(oct)),); name = :x)
dst = LatLonGrid(lat_edges = collect(-90.0:15.0:90.0),
                 lon_edges = collect(0.0:15.0:360.0))

w = conservative_weights(mesh(f), dst)   # build once
out = remap(w, f)                        # DiscreteField{CellLoc} on dst
```

For a real file only the source changes: `fs = load_ugrid("oQU480.ugrid.nc")`
returns a `FieldSet` on an `UnstructuredMesh`, and then
`w = conservative_weights(mesh(fs), dst)` / `out = remap(w, fs[:bottomDepth])`
work as above. Global integrals of the source and the result agree to machine
precision as long as the two meshes cover the same region —
`conservative_weights` verifies that at construction, and the check is
deliberate: an uncovered destination cell is an error, not a silent zero, so a
*global* destination is rejected for a partial-coverage source. That bites for
the tutorial file: `oQU480.ugrid.nc` is ocean-only
(`Σ cell_volume ≈ 0.70 × 4π`), so the global `dst` above fails with
"destination cell 1 overlap sum 0.0 != cell area …". Choose a destination over
the source footprint instead. That file ships with the `ManifoldFields.jl`
tutorial data (`ManifoldFields.jl/examples/data/oQU480.ugrid.nc`), so this part
is prose — the snippet above is self-contained.

**Package status:** private/development package; not registered. See `ROADMAP.md`.
