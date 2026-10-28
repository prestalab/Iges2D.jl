using Iges2D

length(ARGS) == 1 || error("usage: julia --project=. test/cad_export.jl OUTPUT_DIRECTORY")
out = only(ARGS)
mkpath(out)

# All dimensions are millimetres. Non-unit knot domains deliberately test fidelity.
cases = (
    ("line", Line((1, 2), (4, 6))),
    ("quarter_arc", CircularArc((0, 0), (10, 0), (0, 10))),
    ("clockwise_arc", CircularArc((0, 0), (0, 10), (10, 0); ccw=false)),
    ("circle", CircularArc((3, -2), 5.0, 0.0, 2π)),
    ("rational_quarter", NurbsCurve(2, [2, 2, 2, 5, 5, 5],
        [(10, 0), (10, 10), (0, 10)], [1, sqrt(0.5), 1])),
    ("polynomial", NurbsCurve(2, [2, 2, 2, 5, 5, 5],
        [(0, 0), (5, 10), (10, 0)])),
    ("nurbs_circle", as_nurbs(CircularArc((3, -2), 5.0, 0.0, 2π))),
)
for (name, geometry) in cases
    file = IGESFile(; units=:mm, product="Iges2D.jl CAD smoke", date="20261028.120000",
        model_date="20261028.120000")
    add!(file, geometry)
    path = joinpath(out, name * ".igs")
    bytes = write(path, file)
    println(name, ": ", bytes, " bytes")
end

# Exercise metadata continuation with literal punctuation and significant spaces.
metadata = IGESFile(; product=repeat("literal,semicolon; ", 12),
    author="  Iges2D.jl CAD import check  ", organization=repeat("org ", 30))
add!(metadata, last(cases[1]))
write(joinpath(out, "long_metadata.igs"), metadata)
println("Exported ", length(cases) + 1, " IGES files to ", out)
