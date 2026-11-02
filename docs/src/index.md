```@meta
CurrentModule = Iges2D
```

# Iges2D.jl

Write planar CAD curves as ASCII IGES 5.3 files. All geometry lies in the XY
plane (`z = 0`); there is no tessellation of circular arcs or NURBS curves.

| Julia type | IGES entity |
| --- | --- |
| [`CircularArc`](@ref) | 100 — circular arc |
| [`Line`](@ref) | 110 — line segment |
| [`NurbsCurve`](@ref) | 126 — rational B-spline curve |

Requires Julia 1.9 or later. The core writer uses Julia standard libraries;
[BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) integration is
optional and loads automatically when both packages are in the environment.

## Installation

```julia
using Pkg
Pkg.add("Iges2D")
```

## Exporting curves

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
```

See the [repository README](https://github.com/prestalab/Iges2D.jl) for the
arc direction semantics, interpolation, metadata limits and the BSplineKit
integration guide.
