# ASCII IGES 5.3: 72 data columns (64 in P), section letter, seven-digit sequence.
using Dates

const UNITS = Dict(:inch => (1, "IN"), :mm => (2, "MM"), :cm => (10, "CM"),
                   :m => (6, "M"), :um => (9, "MICRON"))

"""
    IGESFile(path=nothing; units=:mm, product="Iges2D.jl", native=product,
             receiver="", preprocessor="Iges2D.jl <version>", author="",
             organization="", date=nothing, model_date=nothing)

A collection of planar curves. Metadata and the exported basename must be
printable ASCII; punctuation and whitespace are preserved in Hollerith strings.
Timestamps are `YYYYMMDD.HHMMSS`; omitted dates use UTC.
"""
mutable struct IGESFile
    path::Union{Nothing,String}
    units::Symbol
    product::String
    native::String
    receiver::String
    preprocessor::String
    author::String
    organization::String
    date::String
    model_date::String
    entities::Vector{Union{Line,CircularArc,NurbsCurve}}
end

function IGESFile(path=nothing; units::Symbol=:mm, product::AbstractString="Iges2D.jl",
                  native::AbstractString=product, receiver::AbstractString="",
                  preprocessor::AbstractString="Iges2D.jl $(pkgversion(@__MODULE__))",
                  author::AbstractString="", organization::AbstractString="",
                  date=nothing, model_date=nothing)
    haskey(UNITS, units) || throw(ArgumentError("unsupported units: $units"))
    timestamp = _timestamp(date)
    IGESFile(isnothing(path) ? nothing : String(path), units, _ascii(product),
             _ascii(native), _ascii(receiver), _ascii(preprocessor), _ascii(author),
             _ascii(organization), timestamp,
             isnothing(model_date) ? timestamp : _timestamp(model_date),
             Union{Line,CircularArc,NurbsCurve}[])
end

"""
    add!(file, entity)

Append a `Line`, `CircularArc` or `NurbsCurve`, returning its one-based entity
index. NURBS knot values and their active parameter interval are preserved.
"""
function add!(f::IGESFile, e)
    e isa Union{Line,CircularArc,NurbsCurve} || throw(ArgumentError("unsupported entity: $(typeof(e))"))
    push!(f.entities, e)
    length(f.entities)
end

"""
    write(io, file::IGESFile)
    write(path, file::IGESFile)
    write(file::IGESFile)

Write CRLF-terminated, 80-byte IGES records and return the byte count. The last
form uses the stored path. The explicit path determines the Global filename.
All records are validated before opening/truncating a path or writing to an IO.
"""
Base.write(io::IO, f::IGESFile) = _write_records(io, _file_records(f, _filename(f)))

function Base.write(path::AbstractString, f::IGESFile)
    records = _file_records(f, basename(path))
    open(path, "w") do io
        _write_records(io, records)
    end
end

function Base.write(f::IGESFile)
    isnothing(f.path) && throw(ArgumentError("no stored path; use write(path, file)"))
    write(f.path, f)
end

function _write_records(io::IO, records)
    bytes = 0
    for record in Iterators.flatten(records)
        bytes += write(io, record, "\r\n")
    end
    bytes
end

function _file_records(f::IGESFile, filename::AbstractString)
    srecs = [_record("IGES 5.3 file written by Iges2D.jl", 'S', 1)]
    grecs = [_record(b, 'G', i) for (i, b) in enumerate(_pack(_global_fields(f, filename); width=72))]
    drecs, precs = String[], String[]
    for (i, e) in enumerate(f.entities)
        de = 2i - 1
        body = _pack(_param_fields(e); width=64)
        offset = length(precs)
        push!(drecs, _record(_dir_line1(e, offset + 1), 'D', de))
        push!(drecs, _record(_dir_line2(e, length(body)), 'D', de + 1))
        for (k, b) in enumerate(body)
            push!(precs, _record(rpad(b, 64) * _integer(de, 8), 'P', offset + k))
        end
    end
    counts = string('S', _integer(length(srecs), 7), 'G', _integer(length(grecs), 7),
                    'D', _integer(length(drecs), 7), 'P', _integer(length(precs), 7))
    # Keep sections separate: concatenating every record would copy the entire file.
    (srecs, grecs, drecs, precs, [_record(counts, 'T', 1)])
end

function _global_fields(f::IGESFile, filename::AbstractString)
    haskey(UNITS, f.units) || throw(ArgumentError("unsupported units: $(f.units)"))
    flag, uname = UNITS[f.units]
    ["1H,", "1H;", _hollerith(f.product), _hollerith(filename),
     _hollerith(f.native), _hollerith(f.preprocessor),
     "32", "38", "6", "308", "15", _hollerith(f.receiver),
     _fmt(1.0), string(flag), _hollerith(uname), "8", "0.016",
     _hollerith(_timestamp(f.date)), "0.0001", _fmt(_max_coord(f.entities)),
     _hollerith(f.author), _hollerith(f.organization),
     "11", "0", _hollerith(_timestamp(f.model_date)), ""]
end

_filename(f::IGESFile) = isnothing(f.path) ? "Iges2D.igs" : basename(f.path)
_hollerith(s::AbstractString) = string(ncodeunits(_ascii(s)), 'H', s)

_points(l::Line) = (l.p1, l.p2)
_points(n::NurbsCurve) = n.points
_coordinate_bound(e::Union{Line,NurbsCurve}) = maximum(max(abs(p.x), abs(p.y)) for p in _points(e))
# A full-circle bound also safely encloses any partial arc.
_coordinate_bound(a::CircularArc) = max(abs(a.center.x), abs(a.center.y)) + _norm(a.start - a.center)
_max_coord(entities) = isempty(entities) ? 0.0 : maximum(_coordinate_bound, entities)

_entity_type(::Line) = 110
_entity_type(::CircularArc) = 100
_entity_type(::NurbsCurve) = 126
_dir_line1(e, pstart) = join(_integer(v, 8) for v in (_entity_type(e), pstart, 0, 1, 0, 0, 0, 0)) * "00000000"
_dir_line2(e, nlines) = join(_integer(v, 8) for v in (_entity_type(e), 0, 0, nlines, 0)) * " "^32

function _param_fields(l::Line)
    fields = String["110"]
    append!(fields, map(_fmt, (l.p1.x, l.p1.y, 0.0, l.p2.x, l.p2.y, 0.0)))
end

function _param_fields(a::CircularArc)
    # Entity 100 only expresses CCW traversal; swap endpoints for a CW locus.
    s, e = a.ccw ? (a.start, a.stop) : (a.stop, a.start)
    fields = String["100"]
    append!(fields, map(_fmt, (0.0, a.center.x, a.center.y, s.x, s.y, e.x, e.y)))
end

function _param_fields(n::NurbsCurve)
    # The public vectors can be mutated after construction. Revalidate at export.
    _check_nurbs(n.degree, n.knots, n.points, n.weights)
    v0, v1 = knot_range(n)
    closed = evaluate(n, v0) == evaluate(n, v1)
    polynomial = all(w -> w == first(n.weights), n.weights)
    fields = string.([126, length(n.points) - 1, n.degree, 1, Int(closed), Int(polynomial), 0])
    append!(fields, map(_fmt, n.knots))
    append!(fields, map(_fmt, n.weights))
    for p in n.points
        append!(fields, (_fmt(p.x), _fmt(p.y), "0."))
    end
    append!(fields, (_fmt(v0), _fmt(v1), "0.", "0.", "1."))
    fields
end

# Numeric tokens stay intact. Long Hollerith strings may span full-width records;
# there must be no inserted padding inside their declared character payload.
function _pack(fields; width::Int)
    lines = String[]
    current = ""
    for (i, field) in enumerate(fields)
        token = String(field) * (i == length(fields) ? ";" : ",")
        if ncodeunits(current) + ncodeunits(token) <= width
            current *= token
        else
            if !isempty(current)
                push!(lines, current)
                current = ""
            end
            if ncodeunits(token) > width
                occursin(r"^[0-9]+H", field) || throw(ArgumentError("numeric field exceeds record width"))
                offset = 1
                while ncodeunits(token) - offset + 1 >= width
                    push!(lines, token[offset:(offset + width - 1)])
                    offset += width
                end
                current = token[offset:end]
            else
                current = token
            end
        end
    end
    isempty(current) || push!(lines, current)
    lines
end

function _integer(value::Integer, width::Int)
    0 <= value < 10^width || throw(ArgumentError("integer $value exceeds $width-column field"))
    lpad(string(value), width)
end

function _record(body::AbstractString, letter::Char, seq::Integer)
    _ascii(body)
    ncodeunits(body) <= 72 || throw(ArgumentError("record body exceeds 72 bytes"))
    seq > 0 || throw(ArgumentError("record sequence must be positive"))
    string(rpad(body, 72), letter, _integer(seq, 7))
end

function _fmt(x::Real)
    value = Float64(x)
    isfinite(value) || throw(ArgumentError("cannot write non-finite number"))
    replace(string(value), 'e' => 'E')
end

function _ascii(s::AbstractString)
    all(b -> 0x20 <= b <= 0x7e, codeunits(s)) ||
        throw(ArgumentError("IGES metadata must be printable ASCII"))
    String(s)
end

function _timestamp(d)
    isnothing(d) && return Dates.format(Dates.now(UTC), dateformat"yyyymmdd.HHMMSS")
    s = String(d)
    occursin(r"^[0-9]{8}\.[0-9]{6}$", s) || throw(ArgumentError("timestamp must be YYYYMMDD.HHMMSS"))
    try
        DateTime(s, dateformat"yyyymmdd.HHMMSS")
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("invalid IGES timestamp: $s"))
    end
    s
end
