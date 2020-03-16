using FASTX
using CodecZlib
using Printf
using BioSequences

mutable struct MerTable{E}
    data::Array{E} # The k-mer stored
    count::Array{Int}   # The number observed
    collision::Int      # Nb. of collision
    unique::Int         # Total number of unique k-mer observed
    capacity::Int       # Size of the arrays
end

function MerTable{E}(capacity) where {E}
    return MerTable(zeros(E, capacity), zeros(Int, capacity),
                    0, 0, capacity)
end

function record!(t::MerTable{E}, m::UInt64) where {E}
    i = m % t.capacity + 1
    t.data[i] = m
    if t.count[i] == 0
        t.unique += 1
    elseif m != t.data[i]
        t.collision += 1
    end
    t.count[i] += 1
end

function Base.show(io::IO, t::MerTable{E}) where {E}
    println(io, "Capacity: ", t.capacity)
    println(io, "Collision: ", t.collision, " (", t.collision / t.capacity * 100, "%)")
    println(io, "Unique (M): ", t.unique / 1e6)
    println(io, "Occupancy: ", t.unique / t.capacity * 100, "%")
end

function build_mertable!(fn, d::MerTable)
    count_seq = 0
    count_k = 0
    for record in open(fn) |> GzipDecompressorStream |> FASTQ.Reader
        s = FASTQ.sequence(record)
        for m in each(DNAMer{31}, s)
            c = convert(UInt64, fwmer(m))
            record!(d, c)
        end
        count_seq += 1
        if count_seq > 1e6
            break
        end
    end
    println(count_k)
    return d
end
