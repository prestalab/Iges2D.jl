module Iges2DBSplineKitExt

import Iges2D: Iges2D, NurbsCurve, Point2
import BSplineKit: Spline, order, knots, basis, coefficients

"""
    NurbsCurve(xs::Spline, ys::Spline)
    NurbsCurve(s::Spline)

Convert two scalar splines sharing a basis, or one spline with two-coordinate
coefficients, to a polynomial NURBS curve. Knots retain their parameter domain.
Only continuous curves (interior multiplicity no greater than degree) are supported.
"""
function Iges2D.NurbsCurve(xs::Spline, ys::Spline)
    bx = basis(xs)
    bx == basis(ys) || throw(ArgumentError("coordinate splines must share the same basis"))
    pts = [Point2(x, y) for (x, y) in zip(coefficients(xs), coefficients(ys))]
    NurbsCurve(order(bx) - 1, collect(knots(bx)), pts)
end

function Iges2D.NurbsCurve(s::Spline)
    coefs = coefficients(s)
    eltype(coefs) <: Number && throw(ArgumentError("scalar spline: pass x and y splines"))
    pts = [Iges2D.point2(c) for c in coefs]
    NurbsCurve(order(basis(s)) - 1, collect(knots(basis(s))), pts)
end

end # module
