# Independent tokenizer for Global fields, including Hollerith data containing
# delimiters, whitespace and physical record boundaries.
function global_values(records)
    text = join(line[1:72] for line in records['G'])
    values = String[]
    i = 1
    while i <= lastindex(text)
        while text[i] == ' '
            i += 1
        end
        j = i
        while isdigit(text[j])
            j += 1
        end
        if j > i && text[j] == 'H'
            n = parse(Int, text[i:j-1])
            push!(values, text[j+1:j+n])
            i = j + n + 1
            while text[i] == ' '
                i += 1
            end
        else
            j = something(findnext(c -> c == ',' || c == ';', text, i))
            push!(values, strip(text[i:j-1]))
            i = j
        end
        @test text[i] == ',' || text[i] == ';'
        text[i] == ';' && break
        i += 1
    end
    values
end

function written_records(file)
    io = IOBuffer()
    count = write(io, file)
    text = String(take!(io))
    @test count == ncodeunits(text)
    iges_records(text)
end

@testset "IGES writer regressions" begin
    empty = written_records(IGESFile(; date="20261028.120000"))
    @test isempty(empty['D']) && isempty(empty['P'])
    @test length(global_values(empty)) == 26

    for (weights, flag) in (([1.0, sqrt(0.5), 1.0], 0), ([2.0, 2.0, 2.0], 1))
        c = NurbsCurve(2, [2, 2, 2, 5, 5, 5], [(10, 0), (10, 10), (0, 10)], weights)
        file = IGESFile()
        add!(file, c)
        records = written_records(file)
        params = only(entity_records(records))[2]
        @test parse(Int, params[6]) == flag
        @test parse.(Float64, params[8:13]) == [2, 2, 2, 5, 5, 5]
        @test parse.(Float64, params[end-4:end]) == [2, 5, 0, 0, 1]
        @test c.knots == [2, 2, 2, 5, 5, 5]
    end

    full = as_nurbs(CircularArc((0, 10), 5, 0, 2π))
    f = IGESFile()
    add!(f, full)
    params = only(entity_records(written_records(f)))[2]
    @test parse(Int, params[5]) == 1                  # closed
    arcfile = IGESFile()
    add!(arcfile, CircularArc((0, 10), 5, 0, 2π))
    records = written_records(arcfile)
    @test parse(Float64, global_values(records)[20]) >= 15
    arcp = only(entity_records(records))[2]
    @test arcp[5:6] == arcp[7:8]                   # exact full-circle endpoints

    # These lengths exercise both sides of physical-record wrapping boundaries.
    for n in (67, 68, 69, 72, 139, 240)
        value = "  ,;" * repeat("q", n) * ";,  "
        file = IGESFile(; product=value, author="  A; B, C  ")
        add!(file, Line((0, 0), (1, 1)))
        globaldata = global_values(written_records(file))
        @test globaldata[3] == value
        @test globaldata[21] == "  A; B, C  "
    end
    @test_throws ArgumentError IGESFile(; product="геометрия")
    @test_throws ArgumentError IGESFile(; author="a\nb")
    @test_throws ArgumentError IGESFile(; date="20260230.120000")
    @test_throws ArgumentError IGESFile(; date="20261028.250000")

    mktempdir() do dir
        file = IGESFile("unused.igs")
        add!(file, Line((0, 0), (1, 1)))
        actual = joinpath(dir, "actual, file.igs")
        @test write(actual, file) == filesize(actual)
        @test global_values(iges_records(read(actual, String)))[4] == basename(actual)
        @test file.path == "unused.igs"

        # Invalid public vector mutation must not clobber an existing destination.
        curve = NurbsCurve(1, [0, 0, 1, 1], [(0, 0), (1, 1)])
        broken = IGESFile()
        add!(broken, curve)
        curve.weights[1] = 0
        original = read(actual)
        @test_throws ArgumentError write(actual, broken)
        @test read(actual) == original
        io = IOBuffer()
        write(io, "prefix")
        @test_throws ArgumentError write(io, broken)
        @test String(take!(io)) == "prefix"
    end
    @test_throws ArgumentError write(IGESFile())
end
