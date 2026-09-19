# Roadmap

`ManifoldRegrid.jl` is a private/development package. Registration is deferred
until the weights API stabilizes.

## Release prerequisites

Tag the siblings **first**, in this order:

1. `ManifoldMeshes` **0.8.0** — `polygon.jl` (clipping, `cell_ring`, area) and
   the id-order batch accessors.
2. `ManifoldFields` **0.1.1** — `location_axis` and `interpolate(f, dest_mesh)`.
3. `ManifoldRegrid` **0.1.0** — last, and only once the two above are tagged.

`ManifoldRegrid` is not consumable until both siblings are released: its
`[compat]` requires `ManifoldMeshes 0.8` (0.7.x has no `cell_ring`, so an
install-and-load against 0.7 dies with `UndefVarError: cell_ring not defined`)
and it calls `ManifoldFields.location_axis`. Until then, `Pkg.develop` on the
sibling paths (or a `rev`-pinned clone) is required.

Accordingly, `.github/workflows/CI.yml`'s "Add sibling packages" step installs
both siblings from GitHub at **pinned revisions** — the merge commits on `main`
that carry the required API — so CI is reproducible and records the exact
sibling generation this package was tested against. The pins are a temporary
bridge: drop them in favour of the released versions once `0.8.0` and `0.1.1`
are tagged.

## Delivered (2026-09)

- First-order (piecewise-constant) conservative remapping:
  `conservative_weights` builds the sparse overlap-area operator from
  spherical Sutherland–Hodgman clipping; `remap(Conservative(), f, dest)`
  applies it with exact destination-area row normalization.
- Bilinear nodal remapping: `remap(Bilinear(), f, dest)`, backed by
  `ManifoldFields.interpolate(f, dest)`.
- Works for every `ManifoldMeshes` grid type, including `UnstructuredMesh`
  sources and destinations; `HEALPixGrid` is the one exception (see Next).

## Next

1. **Second-order reconstruction.** Linear reconstruction inside source cells
   with the gradient term integrated over the same intersection boundary the
   area term uses (one clipping pass serves both), plus a monotonicity limiter.
   The intersection ring returned by `spherical_polygon_intersection` is the
   shared input.
2. **Weight persistence.** Serialize `ConservativeWeights` (`W`, source and
   destination mesh fingerprints) and read ESMF/SCRIP weight files, so a
   weight matrix can be computed once and reused across runs.
3. **General manifolds.** The geometry layer is written in intrinsic terms —
   `side_of_geodesic` is a logarithmic-chart direction test and the area is a
   Stokes boundary integral evaluated in closed form for constant curvature.
   Supporting varying curvature means replacing the closed form with a
   quadrature along the boundary arcs (`exp`/`log`/parallel transport from
   Manifolds.jl) plus row normalization for exact conservation; the interface
   does not change. See `plan/stokes-line-integral-explainer.md` §7.
4. **Partial coverage.** Uncovered destination cells currently have zero row
   sums and are rejected by the conservation check; a future mode should mark
   them explicitly instead.
5. **HEALPix area consistency (upstream).** `Σ cell_volume` over `HEALPixGrid`
   does not equal `4π`: exact at `nside = 1`, **−0.65 %** at `nside = 2`,
   **+0.05 %/+0.11 %/+0.05 %** at `nside = 4/8/16` — non-monotonic in
   resolution, so the geodesic cells of that grid do not tile the sphere
   exactly and the construction-time conservation check rejects `HEALPixGrid`
   as a remapping partner. This lives in the pre-existing
   `cell_nodes`/`cell_volume` construction in `ManifoldMeshes`, not in the
   polygon-geometry layer added here; LatLon, CubedSphere, ReducedGaussian, and
   `UnstructuredMesh` sum to `4π` to machine precision after the
   `spherical_triangle_area` precision fix.
