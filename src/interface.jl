using Logging
"""
The main output object for rdsignal
"""

Base.@kwdef mutable struct WFDBRecord{T} <: AbstractMatrix{T where T<:Real}
    header::Header{<:UnionSpecVector}
    signal::Matrix{T}
    calculated_checksum::Vector{Int64} = Int64[]
    checksum_calculated::Bool = false


    WFDBRecord(h::Header{S},M::AbstractMatrix{Int64}) where {S<:UnionSpecVector} = new{Int64}(h,M)
    WFDBRecord(h::Header{S},M::AbstractMatrix{Float64}) where {S<:UnionSpecVector} = new{Float64}(h,M)

    WFDBRecord(h::Header{S},M::AbstractMatrix{Int64},cs::Vector{Int64}) where {S<:UnionSpecVector} = new{Int64}(h,M,cs,true)
    WFDBRecord(h::Header{S},M::AbstractMatrix{Float64},cs::Vector{Int64}) where {S<:UnionSpecVector} = new{Float64}(h,M,cs,true)
end
Base.size(s::WFDBRecord) = size(s.signal)
Base.length(s::WFDBRecord) = length(s.signal)
Base.getindex(s::WFDBRecord,i::Int) = getindex(s.signal,i)
Base.getindex(s::WFDBRecord,I::Vararg{Int,2}) = s.signal[I...]
Base.IndexStyle(::WFDBRecord) = IndexLinear()
Base.convert(::Type{Vector{T}},s::WFDBRecord) where {T <: Number} = convert(Vector{T}, vec(s.signal))
# Base.Vector{U}(s::WFDBRecord{T}) where {T<:Real,U<:Real} = Vector{U}(s.signal)

export header
header(s::WFDBRecord) = getfield(s,:header)
nsignals(s::WFDBRecord) = header(s) |> nsignals
initial_value(s::WFDBRecord) = header(s) |> initial_value

const CHECKSUM_MODULO = 65536
function validate_checksum(s::WFDBRecord)
    isempty(s.calculated_checksum) && return
    isvalid = Bool[]
    for (a,b) in zip(s.calculated_checksum, checksum(s.header))
        push!(isvalid, mod(b - a, 65536) == 0 )
    end
    if !all(isvalid)
        throw(error("checksums of channels $(findall(!,isvalid)) are invalid"))
    end
    return nothing
end

#pass calls for Header elements are passed through from Signal to
for n in fieldnames(Header)
    expr = :($n(s::WFDBRecord) = s.header |> $n)
    eval(expr)
end

"""
    rdsignal(h::Header,physical::Bool = true)::Tuple{Int64, Matrix{Int32}}

read the samples described in the header's `signal_specs` field.
Header - a parsed header (.hea) file

physical - optional (defaults to true)
        Specifies whether to return signals in physical units in the
        `p_signal` field (True), or digital units in the `d_signal`
        field (False).
"""
function rdrecord(header_path::String,physical::Bool=true)::WFDBRecord
    header = rdheader(header_path)
    rdsignal(header,physical)
end

rdsignal(header::Header) = rdsignal(header::Header, true)::WFDBRecord
function rdsignal(header::Header{T}, physical::Bool;debug::Bool=false)::WFDBRecord where T <: UnionSpecVector

    checksum_calculated = false
    dir = parentdir(header)
    #samples per frame is part of the spec line and samples per signal is part of
    #the Header line
    specline = signalspecline(header)
    if debug
        io = open("/tmp/WFDBLog.txt","w+")
        logger = SimpleLogger(io,Debug)
        with_logger(logger) do
            @info "opening logger"
            samples = _process(dir, signalspecline(header),samples_per_signal(header))
        end
        flush(io)
        close(io)
    else
        samples = _process(dir, signalspecline(header),samples_per_signal(header))
    end

    if T === SingleSpecVector
        _checksum = checksum(samples, header)
        checksum_calculated = true
    end
    if physical
        samples = float(samples) #TODO: change to convert and check out type of checksum...
        dac!(samples, header)
    end

    return checksum_calculated ? WFDBRecord(header,reshape(samples, nsignals(header), :), _checksum) : WFDBRecord(header,reshape(samples, nsignals(header), :))
end

# function _process(dir::String, v::MultiSpecVector,samplespersignal::UInt32)::Vector{Vector{Int64}}
function _process(dir::String, v::MultiSpecVector,samplespersignal::UInt32)::Vector{Int64}
    N = length(v)
    out = Vector{Int64}(undef,N*samplespersignal)
    @debug "out shape:$(size(out))"
    idxs = 0:length(out)-1
    o = ones(Int64,samplespersignal)
    for (i,vi) in enumerate(v)
        out[i:N:length(out)] = _process(dir, vi,samplespersignal)
    end
    return out
    # return reshape(out,length(out))
    # v .|> x-> _process(dir, x,samplespersignal) |> xs->vec(reduce(hcat,xs)')
end

function _process(dir::String, v::SingleSpecVector,samplespersignal::UInt32)::Vector{Int64}
    #TODO: validation step for format. There shouldn't be multi-format in a single file
    fmt = format(v)[1]
    @debug fmt
    n_samples = sum(samples_per_frame(v) * samplespersignal)
    fname = filename(v)[1]
    fname = joinpath(dir, fname)
    extension = get_extension_symbol(fname)

    uniquespf = samples_per_frame(v) |> unique
    uniform = (length(uniquespf) == 1 && uniquespf[1] == 1)
    !uniform && error("non-unity frame sizes not supported")

    extension == :matlab && return read_matlab(fname)

    if extension == :wfdb
        io = open(fname)
        samples = read_binary(io,v, n_samples, fmt)
        @debug "sample type: $(typeof(samples))"
        @debug "length of samples = $(length(samples))"
        close(io)
        return samples
    end
end

function read_matlab(fname::String)::Vector{Int64}
    samples = matread(fname) |> values |> collect
    if length(samples) > 1
        error("more than one matrix in .mat file")
    end
    return convert(Vector{Int64},vec(samples[1]))
end


"""
    wsignal(header::Header, signal::Vector{Int32})

writes a signal file to disk. All of the information required to write the file needs to be specified by the header
"""
function wsignal(header::Header, signal::AbstractVector{Int64})
    fmt = format(header) |> unique
    if length(format) > 1
        throw(error("multi format header writing is not supported"))
    end
    fmt = fmt[1]
    parent = parendir(header)
    if !(isdir(parent))
        throw(error("directory '$(parent)' does not exist"))
    end
    fnames = filename(header)
    uniquefname = unique(fnames)
    if length(uniquefname) > 1
        throw(error("multi-file writers not supported"))
    end
    uniquefname = uniquefname[1]
    full_path = joinpath(parent, uniquefname)
    if isfile(full_path)
        throw(error("file '$(full_path)' already exists"))
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
    function write_binary(io::IO, header::Header, samples::AbstractVector{Int64}, ::WfdbFormat{<:AbstractStorageFormat})

reads a WaveformDB sample file
io - IOBuffer to write to
header - header struct containing the information necessary to *encode* the file
F - the WaveformDB format to encode the output as
"""
function write_binary end
