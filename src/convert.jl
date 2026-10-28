# Conversions between entities.

"""
    clamped_knots(degree, breaks)

Clamped knot vector: first and last break repeated `degree + 1` times, interior
breaks repeated `degree` times.
"""
function clamped_knots(degree::Integer, breaks::AbstractVector{<:Real})
    degree ≥ 1 || throw(ArgumentError("clamped_knots: degree must be >= 1, got $degree"))
    length(breaks) ≥ 2 || throw(ArgumentError("clamped_knots: need at least 2 breaks"))
    all(isfinite, breaks) || throw(ArgumentError("clamped_knots: breaks must be finite"))
    all(i -> breaks[i] > breaks[i-1], 2:length(breaks)) ||
        throw(ArgumentError("clamped_knots: breaks must be strictly increasing"))
    knots = Float64[]
    for (i, v) in enumerate(breaks)
        mult = (i == 1 || i == length(breaks)) ? degree + 1 : degree
        append!(knots, fill(Float64(v), mult))
    end
    knots
end

"""
    as_nurbs(arc::CircularArc; max_piece = π/2)

Exact rational-quadratic NURBS representation of `arc`.  Each piece is a rational
Bézier with weights `1, cos(Δ/2), 1`, so the sweep `Δ` of one piece must stay
below π to keep the middle weight positive.  The arc is split into
`ceil(sweep / max_piece)` pieces joined at the arc's own junction points by
double knots; a full circle becomes four pieces with weights `1, √2/2, ...`.
"""
function as_nurbs(arc::CircularArc; max_piece::Real = π / 2)
    limit = Float64(max_piece)
    isfinite(limit) && 0 < limit < Float64(π) || throw(ArgumentError("max_piece must be finite and lie in (0, pi)"))
    r = _norm(arc.start - arc.center)
    a0 = polar_angle(arc.start, arc.center)
    θ = _arc_sweep(arc)
    ratio = θ / limit
    isfinite(ratio) && ratio < (typemax(Int) - 4) / 2 || throw(ArgumentError("max_piece requires too many pieces"))
    n = max(1, ceil(Int, ratio))
    Δ = θ / n
    w = cos(Δ / 2)
    dir = arc.ccw ? 1.0 : -1.0
    c = arc.center
    points = Point2[]
    weights = Float64[]
    sizehint!(points, 2n + 1)
    sizehint!(weights, 2n + 1)
    for j in 1:n
        a = a0 + dir * (j - 1) * Δ
        b = a + dir * Δ
        m = (a + b) / 2
        if j == 1                                   # shared start point once only
            push!(points, arc.start)
            push!(weights, 1.0)
        end
        push!(points, c + Point2((r / w) * cos(m), (r / w) * sin(m)))
        push!(weights, w)
        push!(points, j == n ? arc.stop : c + Point2(r * cos(b), r * sin(b)))
        push!(weights, 1.0)
    end
    NurbsCurve(2, clamped_knots(2, range(0.0, 1.0; length = n + 1)), points, weights)
end

# --- curve fitting -----------------------------------------------------------------

"""
    fit_nurbs(points::AbstractVector{Point2}; degree=3)

Global B-spline interpolation (Piegl & Tiller, The NURBS Book, Algorithm A9.1)
through `points` (at least `degree + 1`).  Parameters are chord-length, the knot
vector is built by averaging, and the linear system is solved directly (normal
equations square the condition number and destabilise the fit).  Returns a
`NurbsCurve` with `length(points)` control points that passes through every
sample within round-off.
"""
function fit_nurbs(points::AbstractVector{Point2}; degree::Integer = 3)
    n = length(points)
    1 <= degree < n || throw(ArgumentError("fit_nurbs: require 1 <= degree < number of points"))
    # chord-length parameters on [0, 1]
    ub = zeros(n)
    for i in 2:n
        distance = _norm(points[i] - points[i-1])
        distance > 0 || throw(ArgumentError("fit_nurbs: consecutive samples must be distinct"))
        ub[i] = ub[i-1] + distance
        isfinite(ub[i]) && ub[i] > ub[i-1] ||
            throw(ArgumentError("fit_nurbs: chord parameters are not representable at this scale"))
    end
    ub ./= ub[end]
    all(i -> ub[i] > ub[i-1], 2:n) || throw(ArgumentError("fit_nurbs: normalized chord parameters must increase"))
    p = Int(degree)
    # knots by averaging (P&T A9.1): interior knots = mean of p consecutive ubar values
    m = n - 1
    knots = vcat(zeros(p + 1), zeros(m - p), ones(p + 1))
    for j in 1:(m - p)
        knots[p + 1 + j] = sum(ub[(j+1):(j+p)]) / p
    end
    # interpolation matrix (P&T A2.2/A2.3: findSpan + basisFuns, 1-based port)
    A = zeros(n, n)
    for i in 1:n
        t = ub[i]
        # find span: 1-based index of the knot interval containing t
        span = clamp(searchsortedlast(knots, t), p + 1, n)
        N = _basis_funs(span, t, knots, p)
        for j in 1:(p + 1)
            A[i, span - p + j - 1] = N[j, p + 1]
        end
    end
    coordinates = Matrix{Float64}(undef, n, 2)
    for i in eachindex(points)
        coordinates[i, 1] = points[i].x
        coordinates[i, 2] = points[i].y
    end
    controls = A \ coordinates
    NurbsCurve(p, knots, [Point2(controls[k, 1], controls[k, 2]) for k in 1:n], ones(n))
end


"1-based port of P&T A2.3 BasisFuns: N[r+1, j+1] = N_{span-p+r, j}(u)."
function _basis_funs(span::Int, u::Float64, knots::Vector{Float64}, p::Int)
    N = zeros(p + 1, p + 1)
    N[1, 1] = 1.0
    left = zeros(p + 1)
    right = zeros(p + 1)
    for j in 1:p
        left[j+1] = u - knots[span+1-j]
        right[j+1] = knots[span+j] - u
        saved = 0.0
        for r in 0:(j-1)
            temp = N[r+1, j] / (right[r+2] + left[j-r+1])
            N[r+1, j+1] = saved + right[r+2] * temp
            saved = left[j-r+1] * temp
        end
        N[j+1, j+1] = saved
    end
    return N
end
