# Geometry entities (planar) and the checks needed to emit valid IGES parameter
# data.  Everything is placed in the z = 0 model plane.


"""
    Point2(x, y)

Planar point.
"""
struct Point2
    x::Float64
    y::Float64

    function Point2(x::Real, y::Real)
        a, b = Float64(x), Float64(y)
        isfinite(a) && isfinite(b) || throw(ArgumentError("point coordinates must be finite Float64 values"))
        new(a, b)
    end
end

Base.:+(a::Point2, b::Point2) = Point2(a.x + b.x, a.y + b.y)
Base.:-(a::Point2, b::Point2) = Point2(a.x - b.x, a.y - b.y)
Base.:*(s::Real, p::Point2) = Point2(s * p.x, s * p.y)
Base.:*(p::Point2, s::Real) = s * p
Base.:(==)(a::Point2, b::Point2) = a.x == b.x && a.y == b.y
Base.isapprox(a::Point2, b::Point2; kwargs...) =
    isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...)

_norm(p::Point2) = hypot(p.x, p.y)
polar_angle(p::Point2, about::Point2) = atan(p.y - about.y, p.x - about.x)

point2(p::Point2) = p
point2(p::Union{Tuple,AbstractVector}) =
    length(p) == 2 ? Point2(p[1], p[2]) :
    throw(ArgumentError("expected a 2-element point, got $(repr(p))"))
point2(p) = throw(ArgumentError("expected a Point2 or a 2-element tuple/vector, got $(repr(p))"))

"""
    Line(p1, p2)
    Line(x1, y1, x2, y2)

Straight line segment (IGES entity 110).
"""
struct Line
    p1::Point2
    p2::Point2

    Line(p1, p2) = new(point2(p1), point2(p2))
end
Line(x1::Real, y1::Real, x2::Real, y2::Real) = Line(Point2(x1, y1), Point2(x2, y2))

"""
    CircularArc(center, start, stop; ccw = true)
    CircularArc(center, radius, start_angle, stop_angle; ccw = true)

Circular arc in the z = 0 plane (IGES entity 100). Angles are radians. Equal
endpoints or angles differing by whole turns denote a full circle. Small nonzero
arcs remain arcs. `ccw=false` walks clockwise from start to stop; entity 100
preserves the locus by swapping the endpoints on export.
"""
struct CircularArc
    center::Point2
    start::Point2
    stop::Point2
    ccw::Bool

    function CircularArc(c::Point2, p::Point2, q::Point2, ccw::Bool)
        r1, r2 = _norm(p - c), _norm(q - c)
        isfinite(r1) && isfinite(r2) && r1 > 0 && r2 > 0 ||
            throw(ArgumentError("arc radii must be positive and finite"))
        isapprox(r1, r2; rtol=1e-9, atol=0) ||
            throw(ArgumentError("arc endpoints have different radii: $r1 and $r2"))
        new(c, p, q, ccw)
    end
end

CircularArc(center, start, stop; ccw::Bool=true) =
    CircularArc(point2(center), point2(start), point2(stop), ccw)

function CircularArc(center, radius::Real, start_angle::Real, stop_angle::Real;
                     ccw::Bool=true)
    r, a, b = Float64(radius), Float64(start_angle), Float64(stop_angle)
    isfinite(r) && r > 0 || throw(ArgumentError("arc radius must be positive and finite"))
    isfinite(a) && isfinite(b) || throw(ArgumentError("arc angles must be finite"))
    delta = b - a
    isfinite(delta) || throw(ArgumentError("arc angle difference exceeds Float64 range"))
    c = point2(center)
    at(t) = Point2(c.x + r * cos(t), c.y + r * sin(t))
    p = at(a)
    full = iszero(delta) || iszero(rem(delta, 2π))
    q = full ? p : at(b)
    !full && p == q && throw(ArgumentError("distinct arc endpoints are not representable at this scale"))
    CircularArc(c, p, q, ccw)
end

# atan(cross, dot) retains tiny sweeps without subtracting nearly equal angles.
function _arc_sweep(a::CircularArc)
    a.start == a.stop && return 2π
    r = _norm(a.start - a.center)
    ux, uy = (a.start.x - a.center.x) / r, (a.start.y - a.center.y) / r
    vx, vy = (a.stop.x - a.center.x) / r, (a.stop.y - a.center.y) / r
    angle = atan(ux * vy - uy * vx, ux * vx + uy * vy)
    iszero(angle) && throw(ArgumentError("arc sweep is not representable at this scale"))
    mod(a.ccw ? angle : -angle, 2π)
end

"""
    NurbsCurve(degree, knots, points, weights)
    NurbsCurve(degree, knots, points)

Rational B-spline curve (IGES entity 126).  `knots` must hold
`length(points) + degree + 1` non-decreasing values, weights are strictly
positive, and interior knots may not repeat more than `degree` times.  The
parameter range is `[knots[degree+1], knots[end-degree]]`.
"""
struct NurbsCurve
    degree::Int
    knots::Vector{Float64}
    points::Vector{Point2}
    weights::Vector{Float64}

    function NurbsCurve(degree::Integer, knots::Vector{Float64},
                        points::Vector{Point2}, weights::Vector{Float64})
        _check_nurbs(degree, knots, points, weights)
        new(Int(degree), knots, points, weights)
    end
end

NurbsCurve(degree::Integer, knots::AbstractVector{<:Real}, points,
           weights::AbstractVector{<:Real}) =
    NurbsCurve(degree, Float64[k for k in knots], Point2[point2(q) for q in points],
               Float64[w for w in weights])

NurbsCurve(degree::Integer, knots::AbstractVector{<:Real}, points) =
    NurbsCurve(degree, knots, points, ones(Float64, length(points)))

function _check_nurbs(degree::Integer, knots::Vector{Float64}, points::Vector{Point2},
                      weights::Vector{Float64})
    degree ≥ 1 || throw(ArgumentError("nurbs curve: degree must be >= 1, got $degree"))
    n = length(points)
    n > degree ||
        throw(ArgumentError("nurbs curve: degree $degree needs more than $degree control points, got $n"))
    length(weights) == n ||
        throw(ArgumentError("nurbs curve: $n control points but $(length(weights)) weights"))
    length(knots) == n + degree + 1 ||
        throw(ArgumentError("nurbs curve: degree $degree with $n control points needs " *
                            "$(n + degree + 1) knots, got $(length(knots))"))
    all(isfinite, knots) || throw(ArgumentError("nurbs curve: knots must be finite"))
    issorted(knots) || throw(ArgumentError("nurbs curve: knots must be non-decreasing"))
    all(>(0), weights) || throw(ArgumentError("nurbs curve: weights must be strictly positive"))
    all(isfinite, weights) || throw(ArgumentError("nurbs curve: weights must be finite"))
    all(p -> isfinite(p.x) && isfinite(p.y), points) ||
        throw(ArgumentError("nurbs curve: control points must be finite"))
    _check_multiplicities(knots, Int(degree))
    span = knots[n+1] - knots[degree+1]
    isfinite(span) && span > 0 || throw(ArgumentError("nurbs curve: active knot interval must be positive and finite"))
    isfinite(knots[end] - knots[1]) || throw(ArgumentError("nurbs curve: knot span exceeds Float64 range"))
    nothing
end

# This package represents continuous curves: interior multiplicity <= degree.
function _check_multiplicities(knots::Vector{Float64}, degree::Int)
    first = 1
    for i in 2:(length(knots)+1)
        if i > length(knots) || knots[i] != knots[first]
            limit = first == 1 || i > length(knots) ? degree + 1 : degree
            i - first <= limit || throw(ArgumentError("knot $(knots[first]) exceeds multiplicity $limit"))
            first = i
        end
    end
    nothing
end

"""
    normalized(curve::NurbsCurve)

Affinely rescale the active parameter interval to `[0, 1]`. Control points and
weights are unchanged. Normalization is optional; export preserves the domain.
"""
function normalized(n::NurbsCurve)
    a, b = knot_range(n)
    a < b || throw(ArgumentError("nurbs curve: empty parameter range [$a, $b]"))
    (a == 0 && b == 1) && return n
    NurbsCurve(n.degree, (n.knots .- a) ./ (b - a), n.points, n.weights)
end

"Parameter range `[knots[degree+1], knots[end-degree]]` of a curve."
function knot_range(n::NurbsCurve)
    (n.knots[n.degree+1], n.knots[end-n.degree])
end

"""
    evaluate(curve::NurbsCurve, u)

Point of `curve` at parameter `u` (de Boor on homogeneous coordinates).
Parameters outside the curve range are clamped to it.
"""
function evaluate(n::NurbsCurve, u::Real)
    p = n.degree
    isfinite(u) || throw(ArgumentError("nurbs curve: parameter must be finite, got $u"))
    t = Float64(u)
    lo, hi = knot_range(n)
    t = clamp(t, lo, hi)
    k = _find_span(n, t)
    # Only degree+1 control points contribute on this span.
    d = Matrix{Float64}(undef, 3, p + 1)
    weight_scale = maximum(@view n.weights[(k-p):k])
    for j in 0:p
        i = k - p + j
        w = n.weights[i] / weight_scale
        d[1, j+1] = w * n.points[i].x
        d[2, j+1] = w * n.points[i].y
        d[3, j+1] = w
    end
    for r in 1:p, j in p:-1:r
        i = k - p + j
        denom = n.knots[i+p-r+1] - n.knots[i]
        α = (t - n.knots[i]) / denom
        for coordinate in 1:3
            d[coordinate, j+1] = (1 - α) * d[coordinate, j] + α * d[coordinate, j+1]
        end
    end
    Point2(d[1, p+1] / d[3, p+1], d[2, p+1] / d[3, p+1])
end

function _find_span(n::NurbsCurve, t::Float64)
    lo = n.degree + 1
    hi = length(n.points)
    t >= n.knots[hi+1] && return hi
    clamp(searchsortedlast(n.knots, t), lo, hi)
end
