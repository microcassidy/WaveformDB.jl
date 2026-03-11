using WaveformDB: write_binary, parse_signal_spec_line
function setup_writer_tests(fmts)
    D = Dict()

    for fmt in fmts
        ln = lines_mapping[Symbol(fmt)]

        spec_line = Vector{SignalSpecLine}([parse_signal_spec_line(spec_lines[ln])])
        D[fmt] = H(spec_line)
    end

    function writer_test(fmt)
        header = D[fmt]
        fmt = format(header)[1]

        samples = rdsignal(header, false)

        fname = unique(filename(header))
        path = joinpath(parentdir(header), fname[1])
        expectations = read(path)
        io = IOBuffer()
        write_binary(io, samples, fmt)
        seekstart(io)
        reality = read(io)
        close(io)
        @test expectations == reality
    end
    return writer_test
end

@testset "writers" begin
    mt = methods(WaveformDB.read_binary)
    writertest = setup_writer_tests(ty)
    for fmt in ty
        @testset "$fmt" begin
            writertest(fmt)
        end
    end
end
