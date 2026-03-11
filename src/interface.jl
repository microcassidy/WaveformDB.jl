"""
    rdsignal(h::Header,physical::Bool = true)::Tuple{Int64, Matrix{Int32}}

read the samples described in the header's `signal_specs` field.
Header - a parsed header (.hea) file

physical - optional (defaults to true)
        Specifies whether to return signals in physical units in the
        `p_signal` field (True), or digital units in the `d_signal`
        field (False).
"""
struct Signal{T} <: AbstractMatrix{T where T<:Real}
    header::Header{<:UnionSpecVector}
    signal::Matrix{T}
    Signal(h::Header,M::AbstractMatrix{T}) where {T<:Real} = new{T}(h,M)
end

for n in fieldnames(Header)
    expr = :($n(s::Signal) = s.header |> $n)
    eval(expr)
end





Base.size(s::Signal) = size(s.signal)
Base.getindex(s::Signal,i::Int) = s.signal[i]
Base.getindex(s::Signal,I::Vararg{Int,2}) = s.signal[I[1],I[2]]

function rdsignal2(header_path::String,physical::Bool=true)::Signal
    header = rdheader(header_path)
    cs,s = rdsignal(header,physical)
    return Signal(header,s)
end

header(s::Signal) = getfield(s,:header)




rdsignal(header::Header) = rdsignal(header::Header, true)
function rdsignal(header::Header{T}, physical::Bool) where T <: UnionSpecVector
    dir = parentdir(header)

    #samples per frame is part of the spec line and samples per signal is part of
    #the Header line
    specline = signalspecline(header)
    samples = _process(dir, signalspecline(header),samples_per_signal(header))

    if T === SingleSpecVector
        _checksum = checksum(samples, header)
    else
        _checksum = nothing
    end
    if physical
        samples = float(samples) #TODO: change to convert and check out type of checksum...
        dac!(samples, header)
    end
    return _checksum, reshape(samples, nsignals(header), :)
end

function _process(dir::String, v::MultiSpecVector,samplespersignal::UInt32)::Vector{Vector{Int32}}
    v .|> x-> _process(dir, x,samplespersignal)
end

function _process(dir::String, v::SingleSpecVector,samplespersignal::UInt32)::Vector{Int32}
    #TODO: validation step for format. There shouldn't be multi-format in a single file
    fmt = format(v)[1]
    n_samples = sum(samples_per_frame(v) * samplespersignal)
    fname = filename(v)[1]
    fname = joinpath(dir, fname)
    extension = get_extension_symbol(fname)

    uniquespf = samples_per_frame(v) |> unique
    uniform = length(uniquespf) == 1 & uniquespf[1] == 1
    !uniform && error("non-unity frame sizes not supported")

    extension == :matlab && return read_matlab(fname)

    if extension == :wfdb
        io = open(fname)
        samples = read_binary(io,v, n_samples, fmt)
        close(io)
        return samples
    end
end

function read_matlab(fname::String)::Vector{Int32}
    samples = matread(fname) |> values |> collect
    if length(samples) > 1
        error("more than one matrix in .mat file")
    end
    return convert(Vector{Int32},vec(samples[1]))
end


"""
    wsignal(header::Header, signal::Vector{Int32})

writes a signal file to disk. All of the information required to write the file needs to be specified by the header
"""
function wsignal(header::Header, signal::Vector{Int32})
    fmt = format(header) |> unique
    if length(format) > 1
        error("multi format header writing is not supported")
    end
    fmt = fmt[1]
    parent = parendir(header)
    if !(isdir(parent))
        error("directory '$(parent)' does not exist")
    end
    fnames = filename(header)
    uniquefname = unique(fnames)
    if length(uniquefname) > 1
        error("multi-file writers not supported")
    end
    uniquefname = uniquefname[1]
    full_path = joinpath(parent, uniquefname)
    if isfile(full_path)
        error("file '$(full_path)' already exists")
    end
    open(full_path, 'w') do io
        write_binary(io, header, samples, fmt)
    end
end

"""
    read_binary(io::IO, header::Header, ::WfdbFormat{<:AbstractStorargeFormat})

reads a WaveformDB sample file
io - IOBuffer for the file ("eg 100.dat")
header - header struct containing the information necessary to *decode* the file
F - the format that the file is in
"""
function read_binary end

"""
    function write_binary(io::IO, header::Header, samples::Vector{Int32}, ::WfdbFormat{<:AbstractStorageFormat})

reads a WaveformDB sample file
io - IOBuffer to write to
header - header struct containing the information necessary to *encode* the file
F - the WaveformDB format to encode the output as
"""
function write_binary end
