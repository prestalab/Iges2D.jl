using Iges2D
using Test

Bernstein(u, i, n) = binomial(n, i) * u^i * (1 - u)^(n - i)

# --- independent structural checker (written from the IGES 5.3 spec) --------

function iges_records(text)
    @test isascii(text)
    lines = split(text, "\r\n"; keepempty = true)
    @test pop!(lines) == ""
    @test all(l -> ncodeunits(l) == 80, lines)
    secs = Dict(c => String[] for c in ('S', 'G', 'D', 'P', 'T'))
    for l in lines
        push!(secs[l[73]], l)
    end
    @test [l[73] for l in lines] == [c for c in ('S', 'G', 'D', 'P', 'T') for _ in secs[c]]
    check_structure(secs)
    secs
end

function check_structure(secs)
    @test length(secs['T']) == 1
    t = only(secs['T'])
    for (i, c) in enumerate(('S', 'G', 'D', 'P'))
        @test t[8(i-1)+1] == c
        @test length(secs[c]) == parse(Int, strip(t[8(i-1)+2:8i]))
    end
end

function _sequence_ok(secs)
    for (c, recs) in secs
        all(((i, r),) -> parse(Int, strip(r[74:80])) == i, enumerate(recs)) ||
            return false
    end
    true
end

function entity_records(secs)
    out = Vector{Tuple{Int,Vector{String}}}()   # (entity type, parameter strings)
    d, p = secs['D'], secs['P']
    @test length(d) % 2 == 0
    for i in 1:2:length(d)
        first, second = d[i], d[i+1]
        etype = parse(Int, strip(first[1:8]))
        pstart = parse(Int, strip(first[9:16]))
        nlines = parse(Int, strip(second[25:32]))
        de = parse(Int, strip(first[74:80]))
        @test parse(Int, strip(second[74:80])) == de + 1
        @test 1 <= pstart <= length(p)
        @test pstart + nlines - 1 <= length(p)
        for k in 0:(nlines-1)
            rec = p[pstart+k]
            @test strip(rec[65:72]) == string(de)   # back pointer to the D entry
        end
        body = join(strip(rec[1:64]) for rec in p[pstart:pstart+nlines-1])
        @test endswith(body, ";")
        params = split(chop(body), ',')
        @test parse(Int, params[1]) == etype
        push!(out, (etype, params))
    end
    out
end


@testset "Iges2D.jl" begin
    @testset "constructors and validation" begin
        a = CircularArc((0, 0), 2, 0, π / 2)
        @test a.start ≈ Point2(2.0, 0.0)
        @test a.stop ≈ Point2(0.0, 2.0) atol = 1e-12
        @test CircularArc((0, 0), 2, 0, 2π).start ≈ CircularArc((0, 0), 2, 0, 2π).stop atol = 1e-12
        @test CircularArc((0, 0), 2, 2π, 0; ccw = false).ccw == false
        @test_throws ArgumentError CircularArc((0, 0), (1, 0), (2, 0))
        @test_throws ArgumentError CircularArc((0, 0), 0.0, 0, 1)
        @test_throws ArgumentError CircularArc((0, 0), -1.0, 0, 1)
        @test_throws ArgumentError CircularArc((0, 0), Inf, 0, 1)
        @test_throws ArgumentError Line((0, 0), (1,))
        n = NurbsCurve(2, clamped_knots(2, [0.0, 1.0]), [(0, 0), (1, 1), (2, 0)])
        @test n.weights == [1.0, 1.0, 1.0]
        @test_throws ArgumentError NurbsCurve(2, [0, 0, 0, 1, 1], [(0, 0), (1, 1), (2, 0)], ones(3))
        @test_throws ArgumentError NurbsCurve(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)], [1, 1, 0])
        @test_throws ArgumentError NurbsCurve(2, [0, 0, 1, 0, 1, 1], [(0, 0), (1, 1), (2, 0)], ones(3))
        @test_throws ArgumentError NurbsCurve(1, [0, 0, 1, 1], [(0, 0), (1, 1), (2, 0)], ones(3))
        # Six control points and nine knots: rejection is multiplicity, not count.
        @test_throws ArgumentError NurbsCurve(2, [0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1],
                                              [(i, 0) for i in 0:5], ones(6))
    end

    @testset "evaluate" begin
        # cubic Bézier checked against the Bernstein form
        p = [Point2(0, 0), Point2(1, 2), Point2(3, 2), Point2(4, 0)]
        c = NurbsCurve(3, [0, 0, 0, 0, 1, 1, 1, 1], p)
        for u in (0.0, 0.1, 0.5, 0.9, 1.0)
            b = sum(Bernstein(u, i, 3) * p[i+1] for i in 0:3)
            @test evaluate(c, u) ≈ b
        end
        @test evaluate(c, -1.0) == evaluate(c, 0.0)
        @test evaluate(c, 2.0) == evaluate(c, 1.0)
        # Rational quadratic midpoint: numerator w/2, denominator (1+w)/2.
        q = NurbsCurve(2, [0, 0, 0, 1, 1, 1], [(0, 0), (1, 1), (2, 0)], [1, cos(π / 4), 1])
        m = evaluate(q, 0.5)
        @test m ≈ Point2(1.0, sqrt(2) - 1.0)
    end

    @testset "as_nurbs" begin
        a = CircularArc((1, 2), 5, 0, π / 2)
        n = as_nurbs(a)
        @test n.degree == 2
        @test length(n.points) == 3
        @test n.weights ≈ [1.0, cos(π / 4), 1.0]
        for u in range(0, 1, length = 25)
            p = evaluate(n, u)
            @test hypot((p - a.center).x, (p - a.center).y) ≈ 5.0 atol = 1e-12
        end
        @test evaluate(n, 0.0) ≈ a.start
        @test evaluate(n, 1.0) ≈ a.stop
        # Clockwise traversal keeps the endpoints and takes the complementary arc.
        cw = as_nurbs(CircularArc((1, 2), 5, 0, π / 2; ccw = false))
        @test evaluate(cw, 0.0) ≈ a.start      # CW arc keeps user start/stop
        @test evaluate(cw, 1.0) ≈ a.stop
        @test hypot((evaluate(cw, 0.5) - a.center).x, (evaluate(cw, 0.5) - a.center).y) ≈ 5.0 atol = 1e-12
        # full circle: four quadratic pieces, all points on the circle
        full = as_nurbs(CircularArc((0, 0), 3, 0, 2π))
        @test length(full.points) == 9
        @test full.weights ≈ [1, sqrt(2) / 2, 1, sqrt(2) / 2, 1, sqrt(2) / 2, 1, sqrt(2) / 2, 1]
        for u in range(0, 1, length = 65)
            @test hypot(evaluate(full, u).x, evaluate(full, u).y) ≈ 3.0 atol = 1e-12
        end
        @test evaluate(full, 0.0) ≈ evaluate(full, 1.0) atol = 1e-12
        # obtuse piece angles stay exact through several pieces
        tri = as_nurbs(CircularArc((0, 0), 1, 0, 2.4); max_piece = 0.6)
        @test length(tri.points) == 9
        for u in range(0, 1, length = 41)
            @test hypot(evaluate(tri, u).x, evaluate(tri, u).y) ≈ 1.0 atol = 1e-12
        end
    end

    @testset "file structure" begin
        f = IGESFile("demo.igs"; date = "20240101.120000", model_date = "20240101.120000")
        @test_throws ArgumentError IGESFile(units = :parsec)
        @test_throws ArgumentError IGESFile(date = "1.1.2024")
        @test add!(f, CircularArc((0, 0), 10, 0, π / 2)) == 1
        @test add!(f, Line((0, 0), (3, 4))) == 2
        @test add!(f, as_nurbs(CircularArc((1, 2), 5, 0, π / 2))) == 3
        @test_throws ArgumentError add!(f, (0, 0))
        io = IOBuffer()
        n = write(io, f)
        text = String(take!(io))
        @test n == sizeof(text)
        secs = iges_records(text)
        @test _sequence_ok(secs)
        ents = entity_records(secs)
        @test first.(ents) == [100, 110, 126]

        # entity 100: planar, start/stop order
        (_, params) = ents[1]
        @test length(params) == 8
        @test parse(Float64, params[2]) == 0.0                 # z displacement
        @test parse(Float64, params[3]) ≈ 0.0 && parse(Float64, params[4]) ≈ 0.0
        @test parse(Float64, params[5]) ≈ 10.0 && parse(Float64, params[6]) ≈ 0.0
        @test isapprox(parse(Float64, params[7]), 0.0; atol = 1e-12) && parse(Float64, params[8]) ≈ 10.0

        # entity 110
        (_, params) = ents[2]
        @test length(params) == 7
        @test parse(Float64, params[2]) ≈ 0.0 && parse(Float64, params[3]) ≈ 0.0
        @test parse(Float64, params[5]) ≈ 3.0 && parse(Float64, params[6]) ≈ 4.0

        # entity 126: field layout K, M, flags, knots, weights, points, V, normal
        (_, params) = ents[3]
        n = as_nurbs(CircularArc((1, 2), 5, 0, π / 2))
        K = length(n.points) - 1
        @test parse(Int, params[2]) == K
        @test parse(Int, params[3]) == 2
        @test parse(Int, params[4]) == 1                       # planar
        @test parse(Int, params[5]) == 0                       # open
        @test parse(Int, params[6]) == 0                       # rational (1 = polynomial)
        @test parse(Int, params[7]) == 0                       # non-periodic
        i = 8
        nknots = parse(Int, params[2]) + parse(Int, params[3]) + 2   # K + M + 2
        knots = [parse(Float64, params[j]) for j in i:(i + nknots - 1)]
        i += nknots
        weights = [parse(Float64, params[j]) for j in i:(i + K)]     # K + 1 weights
        i += K + 1
        pts = [Point2(parse(Float64, params[j]), parse(Float64, params[j+1]))
               for j in i:3:(i + 3 * (K + 1) - 1)]
        i += 3 * (K + 1)
        tail = params[i:end]
        @test length(tail) == 5                                # V(0), V(1), normal
        @test parse(Float64, tail[1]) == 0.0                   # V(0)
        @test parse(Float64, tail[2]) == 1.0                   # V(1)
        @test parse(Float64, tail[3]) == 0.0                   # normal (0, 0, 1)
        @test parse(Float64, tail[4]) == 0.0
        @test parse(Float64, tail[5]) == 1.0
        @test knots == n.knots
        @test weights == n.weights
        @test pts == n.points

        # clockwise arcs are swapped to the CCW form entity 100 requires
        f2 = IGESFile(; date = "20240101.120000")
        add!(f2, CircularArc((0, 0), 10, 0, π / 2; ccw = false))
        io2 = IOBuffer()
        write(io2, f2)
        ents2 = entity_records(iges_records(String(take!(io2))))
        (_, params) = ents2[1]
        @test parse(Float64, params[5]) ≈ 0.0 atol = 1e-12   # CCW: (0,10) -> (10,0)
        @test parse(Float64, params[6]) ≈ 10.0
        @test parse(Float64, params[7]) ≈ 10.0
        @test parse(Float64, params[8]) ≈ 0.0 atol = 1e-12
    end


    @testset "BSplineKit extension" begin
        using BSplineKit: BSplineBasis, Spline, BSplineKit

        basis = BSplineBasis(3, [0.0, 0.25, 0.5, 1.0])
        ts = range(0, 1, length = 5)
        xs = Spline(basis, cospi.(2 .* ts))
        ys = Spline(basis, sinpi.(2 .* ts))
        c = NurbsCurve(xs, ys)
        @test c.degree == 2
        for u in (0.0, 0.1, 0.5, 0.9, 1.0)
            @test evaluate(c, u) ≈ Point2(xs(u), ys(u))
        end
        # point-valued coefficients
        pts = Spline(basis, [Point2(cospi(2t), sinpi(2t)) for t in ts])
        c2 = NurbsCurve(pts)
        for u in (0.0, 0.3, 0.7, 1.0)
            @test evaluate(c2, u) ≈ evaluate(c, u)
        end
        other = Spline(BSplineBasis(3, [0.0, 0.5, 1.0]), [1.0, 2.0, 3.0, 4.0])
        @test_throws ArgumentError NurbsCurve(xs, other)
        # unnormalized knot range is rescaled by normalized()
        off = BSplineBasis(3, [2.0, 2.25, 2.5, 3.0])
        zs = Spline(off, [3.0, 1.0, 2.0, 4.0, 5.0])
        ws = Spline(off, [1.0, 4.0, 2.0, 5.0, 3.0])
        c3 = Iges2D.normalized(NurbsCurve(zs, ws))
        @test c3.knots[1] == 0.0 && c3.knots[end] == 1.0
        for u in (0.0, 0.4, 1.0)
            @test evaluate(c3, u) ≈ Point2(zs(2.0 + u), ws(2.0 + u))
        end
    end
end

@testset "fit_nurbs" begin
    # global cubic interpolation through a quarter-circle sample
    pts = [Point2(cos(t), sin(t)) for t in range(0, π / 2; length=9)]
    c = fit_nurbs(pts)
    @test c.degree == 3
    @test length(c.points) == 9
    @test evaluate(c, 0.0) == pts[1]
    @test evaluate(c, 1.0) == pts[end]
    # passes through the samples exactly at their chord-length parameters
    ub = zeros(length(pts))
    for i in 2:length(pts)
        ub[i] = ub[i-1] + hypot(pts[i].x - pts[i-1].x, pts[i].y - pts[i-1].y)
    end
    ub ./= ub[end]
    for i in eachindex(pts)
        @test isapprox(evaluate(c, ub[i]).x, pts[i].x; atol = 5e-4)
        @test isapprox(evaluate(c, ub[i]).y, pts[i].y; atol = 5e-4)
    end
    # close to the true circle between samples
    rdev = maximum(abs(hypot(evaluate(c, u).x, evaluate(c, u).y) - 1.0)
                   for u in range(0, 1; length=200))
    @test rdev < 1e-3
    @test_throws ArgumentError fit_nurbs(pts[1:2])
    @test_throws ArgumentError fit_nurbs(pts; degree=0)
end

include("geometry_regression.jl")
include("writer_regression.jl")
include("quality.jl")
