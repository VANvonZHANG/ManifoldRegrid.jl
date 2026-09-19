# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Package skeleton with the `AbstractRemapScheme` / `Conservative` / `Bilinear`
  scheme tags.
- `conservative_weights` / `ConservativeWeights`: first-order overlap-area
  weights from spherical Sutherland–Hodgman clipping, with k-d tree candidate
  search and construction-time conservation invariants.
- `remap`: `remap(Conservative(), f, dest)`, `remap(Bilinear(), f, dest)`, and
  `remap(w::ConservativeWeights, f)`. Scheme and field location are checked by
  dispatch.
