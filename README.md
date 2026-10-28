# Iges2D.jl

Write planar CAD curves as ASCII IGES 5.3 files. All geometry lies in the XY
plane (`z = 0`); there is no tessellation of circular arcs or NURBS curves.

| Julia type | IGES entity |
| --- | --- |
| `CircularArc` | 100 — circular arc |
| `Line` | 110 — line segment |
| `NurbsCurve` | 126 — rational B-spline curve |

Requires Julia 1.9 or later. The core writer uses Julia standard libraries;
[BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) integration is optional.

## Install

Before General registration is merged, install from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/prestalab/Iges2D.jl")
```

## Export curves

```julia
using Iges2D

file = IGESFile("detail.igs"; units=:mm, product="Example")
add!(file, Line((0, 0), (10, 0)))
add!(file, CircularArc((5, 0), 5.0, 0.0, π))  # centre, radius, angles in radians
add!(file, CircularArc((5, 5), (10, 5), (0, 5)))  # centre, start, stop

# Exact rational quadratic quarter circle.
curve = NurbsCurve(2, [0, 0, 0, 1, 1, 1],
    [(10, 0), (10, 10), (0, 10)], [1, sqrt(0.5), 1])
add!(file, curve)
add!(file, as_nurbs(CircularArc((0, 0), 3.0, 0.0, 2π)))
write(file)

# Explicit destinations are also supported; write returns the byte count.
write("copy.igs", file)
io = IOBuffer()
write(io, file)
```

`add!` returns the entity's one-based index. `Point2`, two-element tuples, and
two-element vectors are accepted as points. Coordinates, knots, and weights are
stored as `Float64`.

### Arc direction and circles

Arcs run counter-clockwise by default. `ccw=false` requests clockwise traversal
from the supplied start to stop. Entity 100 can only express counter-clockwise
traversal: its endpoints are swapped on export, preserving the geometric locus
but reversing the parameter direction. Use `as_nurbs(arc)` to retain traversal
direction in the exported entity 126.

Equal endpoints denote a full circle. The angle constructor also accepts equal
angles or one full turn. `as_nurbs` splits an arc into exact rational quadratic
pieces (default maximum sweep `π/2`). Its parameter is **not** linear in the
polar angle. Interpolation through sampled circle points is not an exact circle.

### NURBS and interpolation

A degree-`p` curve with `n` control points requires `n + p + 1` knots and `n`
positive weights. Omitting the weights gives unit weights. Knots must be finite
and nondecreasing, with a nonempty active parameter interval
`[knots[p+1], knots[end-p]]`. This package supports continuous curves: interior
knot multiplicity is at most `p`.

`evaluate(curve, u)` evaluates in homogeneous coordinates and clamps `u` to the
active interval. Export preserves knot values and the active interval; IGES
does **not** require normalization to `[0, 1]`. Use `normalized(curve)` explicitly
when an affine parameter change is wanted.

```julia
using Iges2D

samples = [Point2(t, sin(t)) for t in range(0, π; length=9)]
interpolant = fit_nurbs(samples; degree=3)
midpoint = evaluate(interpolant, 0.5)
```

`fit_nurbs` performs global polynomial B-spline interpolation with chord-length
parameters and averaged knots. It is not a least-squares fit; consecutive
coincident samples are rejected. `clamped_knots(p, breaks)` repeats endpoint
breaks `p+1` times and interior breaks `p` times, for piecewise Bézier segments;
it is not a simple-knot interpolation helper.

## BSplineKit integration

Install and load BSplineKit to activate the extension:

```julia
using Iges2D
using BSplineKit: BSplineBasis, Spline

basis = BSplineBasis(3, [0.0, 0.25, 0.5, 1.0])  # order 3 means degree 2
xs = Spline(basis, [0.0, 1.0, 2.0, 3.0, 4.0])
ys = Spline(basis, [0.0, 2.0, 1.0, 2.0, 0.0])
curve = NurbsCurve(xs, ys)
```

Coordinate splines must share the same basis. `NurbsCurve(spline)` also accepts
point-valued coefficients. Conversion preserves the knot domain and uses unit
weights; create `NurbsCurve` with explicit weights for rational curves.

## File metadata and limits

- Units: `:mm`, `:inch`, `:cm`, `:m`, `:um`. Coordinates are interpreted in those
  units, not rescaled by the writer.
- Records contain 80 ASCII bytes followed by CRLF. Metadata and the exported
  basename must be printable ASCII. Hollerith strings retain punctuation and spaces.
- `date` and `model_date` accept `"YYYYMMDD.HHMMSS"`; omitted dates use UTC.
- No IGES reader, surfaces, solids, transforms, colours or topology entities.
- Accuracy depends on coordinate scale and floating-point conditioning; no fixed
  absolute geometric tolerance is guaranteed for all inputs.

## Verification

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests cover numerical evaluation, interpolation and IGES record semantics. CI
runs Julia tests on Linux, Windows and macOS, including the minimum Julia version,
and keeps coverage reports. A separate OpenCASCADE import check validates exported
curves without using Iges2D.jl to read them:

```sh
julia --project=. test/cad_export.jl cad-smoke
python -m pip install cadquery-ocp==7.9.3.1.1
python test/cad_import.py cad-smoke
```

Python/OpenCASCADE are verification tools, not dependencies of the Julia package.

## License and provenance

MIT; see [LICENSE](LICENSE). IGES field definitions were cross-checked against
[libIGES](https://github.com/ramcdona/libIGES),
[forIGES](https://github.com/rweed/forIGES) and OpenCASCADE.

Substantial implementation, tests and documentation were developed with AI
assistance. Automated review is not a substitute for human maintainer review;
General registration requires the maintainer to review and understand the code.
