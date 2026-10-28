using BSplineKit: BSplineKit, BSplineBasis, Spline

@testset "Geometry regressions" begin
    @test_throws ArgumentError Point2(Inf, 0)
    @test_throws ArgumentError Point2(0, NaN)
    @test_throws ArgumentError Line((0, 0), (Inf, 0))
    @test_throws ArgumentError CircularArc(Point2(0, 0), Point2(1, 0), Point2(0, 2), true)
    @test_throws ArgumentError CircularArc((Inf, 0), 1, 0, 1)
    @test Line(Point2(0, 0), (1, 2)).p2 == Point2(1, 2)

    circle = CircularArc((3, -2), 5.0, 0.0, 2π)
    @test circle.start == circle.stop
    @test CircularArc((0, 0), 1, 1e300, 1e300; ccw=false).start ==
          CircularArc((0, 0), 1, 1e300, 1e300).stop
    full = as_nurbs(circle)
    @test first(full.points) == last(full.points)
    @test evaluate(full, 0) == evaluate(full, 1)

    small = CircularArc((0, 0), 1, 0, 1e-13)
    short = as_nurbs(small)
    @test short.points[1] != short.points[end]
    @test evaluate(short, 0) == small.start
    @test evaluate(short, 1) == small.stop
    middle = evaluate(short, 0.5)
    @test 0 < middle.y < 1e-13
    @test middle.x ≈ 1
    for invalid in (0.0, -1.0, π, 2π, Inf, NaN)
        @test_throws ArgumentError as_nurbs(small; max_piece=invalid)
    end
    @test_throws ArgumentError clamped_knots(2, [0, 0, 1])
    @test_throws ArgumentError clamped_knots(2, [0, Inf])
    @test_throws ArgumentError NurbsCurve(2, [-2, -1, 0, 0, 1, 2], [(0, 0), (1, 2), (3, 0)])

    # Homogeneous weights can be rescaled without changing the curve.
    weighted = NurbsCurve(2, [0, 0, 0, 1, 1, 1], [(10, 0), (10, 10), (0, 10)],
        [1e308, sqrt(0.5) * 1e308, 1e308])
    @test evaluate(weighted, 0.5) ≈ Point2(sqrt(50), sqrt(50))
    @test_throws ArgumentError evaluate(weighted, NaN)

    # Nonuniform rational curve checked against independent polynomial splines.
    basis = BSplineBasis(4, [-2.0, -0.1, 0.8, 3.0])
    n = length(basis)
    controls = [Point2(sin(i), cos(2i)) for i in 1:n]
    weights = [0.4 + i / n for i in 1:n]
    curve = NurbsCurve(3, collect(BSplineKit.knots(basis)), controls, weights)
    sx = Spline(basis, [controls[i].x * weights[i] for i in 1:n])
    sy = Spline(basis, [controls[i].y * weights[i] for i in 1:n])
    sw = Spline(basis, weights)
    scaled = normalized(curve)
    for t in range(-2, 3; length=31)
        expected = Point2(sx(t) / sw(t), sy(t) / sw(t))
        @test evaluate(curve, t) ≈ expected atol=1e-13
        @test evaluate(scaled, (t + 2) / 5) ≈ expected atol=1e-13
    end

    # Irregular sample spacing exposes span-search and chord-parameter mistakes.
    samples = [Point2(t, sin(t)) for t in (0.0, 0.01, 0.7, 1.1, 2.0, 3.0, 4.0)]
    interpolant = fit_nurbs(samples)
    params = zeros(length(samples))
    for i in 2:length(samples)
        params[i] = params[i-1] + hypot(samples[i].x - samples[i-1].x, samples[i].y - samples[i-1].y)
    end
    params ./= last(params)
    for i in eachindex(samples)
        @test evaluate(interpolant, params[i]) ≈ samples[i] atol=1e-12
    end
    duplicate = copy(samples)
    duplicate[3] = duplicate[2]
    @test_throws ArgumentError fit_nurbs(duplicate)
end
