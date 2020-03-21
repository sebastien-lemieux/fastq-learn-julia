# using FASTX
# using CodecZlib
# using BioSequences

ktype = UInt64    # Up to 32-mer
ctype = UInt32  # Up to 4G occurrences

mutable struct MerTable
    key::Array{ktype}
    count::Array{ctype}
    overflow::Array{ktype}
    n_kc::ctype  # nb. of keys in key-count table
    n_of::ctype  # nb. of keys in overflow
end

function MerTable(overflowsize)
    return MerTable(Array{ktype}(undef, 0), Array{ktype}(undef, 0), Array{ktype}(undef, UInt(overflowsize)),
                    0, 0)
end

m = MerTable(10)

function find_index(t::MerTable, k::ktype, low::Int=1, high::Int=length(t.key))
    # Assumes the key exists
    println((low, high))
    if low == high
        return low
    end
    mid = div(low + high, 2)
    if k <= t.key[mid]
        return find_index(t, k, low, mid)
    else
        return find_index(t, k, mid + 1, high)
    end
end

# m.key = 1:7
# find_index(m, ktype(0))

function merge_of(t::MerTable)
    unique_of = unique(sort(t.overflow))
    new_key = vcat(t.key, unique_of)

end

function record!(t::MerTable, m::ktype)
    # Search in kc
    t.n_of += 1
    t.overflow[t.n_of] = m
    if t.n_of == length(t.overflow)
        merge_of(t)
    end
    return (t.n_kc, t.n_of)
end

function build_mertable!(fn, d::MerTable{E}) where {E}
    count_seq = 0
    count_k = 0
    for record in open(fn) |> GzipDecompressorStream |> FASTQ.Reader
        s = FASTQ.sequence(record)
        for m in each(DNAMer{31}, s)
            c = convert(UInt64, fwmer(m))
            record!(d, c)
        end
        count_seq += 1
        if count_seq > (8e6 / 56)
            break
        end
    end
end
