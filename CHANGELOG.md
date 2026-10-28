# Changelog

## 0.1.0

Initial public release as **Iges2D** (the unregistered development module was
named IGES).

- Export planar circular arcs (100), line segments (110), and NURBS curves (126).
- Correct rational/polynomial and closed-curve flags, termination counts,
  byte widths, and long Hollerith strings; preserve metadata punctuation.
- Preserve NURBS parameter domains instead of silently normalizing on write.
- Validate output before truncating a destination and use its actual filename
  in the Global section.
- Keep full-circle endpoints exactly equal; retain small nonzero arcs and avoid
  angle-adjustment loops. Reject nonfinite geometry and invalid piece sweeps.
- Restrict the de Boor workspace to the contributing degree+1 control points.
- Retain chord-length interpolation with one factorization for both coordinates;
  reject consecutive coincident samples.
- Integrate BSplineKit through optional package extensions and public accessors.
- Add MIT licensing, compatibility bounds, cross-platform CI, Aqua checks,
  coverage artifacts, and an independent OpenCASCADE import check.

Local release checks: full suite on Julia 1.9.4 and 1.12.6, 877 checks on each;
8 OpenCASCADE 7.9.3 import scenarios covering geometry, metadata and parameter
domains. These are tested cases, not a guarantee for every CAD application.
