using CodecZlib

cd("C:\\Users\\slemi\\prog\\Load")
fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"

struct FastqFile
    io::IO
end

function Base.iterate(f::FastqFile, state=1)
    if eof(f.io)
        return nothing
    end

    readline(f.io) # title
    seq = readline(f.io)
    readline(f.io) # +
    readline(f.io) # qual
    return (seq, state+1)
end

charToUint = Dict{Char,UInt64}('A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 0)

mutable struct MerGenerator
    k::UInt64
    mask::UInt64
    seq::String
    l::Int
end

function MerGenerator(k::UInt64, seq::String)
    mask = (UInt64(1) << (k << UInt64(1))) - UInt64(1)
    return MerGenerator(k, mask, seq, length(seq))
 end

@inline function Base.iterate(g::MerGenerator, state=(Int(1),UInt64(0)))
    i, cur = state
    if i > g.l
        return nothing
    end
    while i < g.k
        c = charToUint[g.seq[i]]
        # println("pre",c, cur)
        cur = ((cur << UInt64(2)) + c) & g.mask
        i += 1
    end
    c = charToUint[g.seq[i]]
    # println("post",c)
    cur = ((cur << UInt64(2)) + c) & g.mask
    i += 1
    return cur, (i, cur)
end

Base.eltype(::Type{MerGenerator}) = UInt64

# m = MerGenerator(UInt(5), "AAAAAAAACAA")
# for s in m
#     println(s)
# end

# 566k reads per sec

function count_seq(fn)
    f = FastqFile(open(fn) |> GzipDecompressorStream)

    sumofk::UInt = 0
    count = 0
    for s in f
        count += 1
        m = MerGenerator(UInt(31), s)
        for kmer in m
            sumofk += kmer
        end
        if count == 100000
            break
        end
    end

    println(count)
end

@time count_seq(fn)
