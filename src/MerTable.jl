using FASTX
using CodecZlib
using Printf
using BioSequences

include("common.jl")

mutable struct MerTable
    data::Array{ktype}  # The k-mer stored
    count::Array{ctype} # The number observed
    # stats
    collision::Int      # Nb. of collision
    unique::Int         # Total number of unique k-mer observed
    capacity::Int       # Size of the arrays
end

function MerTable(capacity)
    c = Int(capacity)
    # This init will only work for k <= 31, at 32, the typemax will exist, this
    # will result in an incorrect (-1) count of the uniques.
    return MerTable(fill(typemax(ktype), c), zeros(ctype, c), 0, 0, c)
end

function record!(t::MerTable, m::ktype)
    index = m % t.capacity + 1
    if t.data[index] == m
        t.count[index] += 1
        return
    end
    while t.count[index] > 0
        if t.data[index] == m
            t.count[index] += 1
            return
        end
        t.collision += 1
        index += 1
        if index > t.capacity
            index -= t.capacity
        end
    end
    # Now t.count[index] == 0
    t.data[index] = m
    t.count[index] = 1
    t.unique += 1
end

function Base.show(io::IO, t::MerTable)
    println(io, "Capacity: ", t.capacity)
    println(io, "Collision: ", t.collision, " (", t.collision / t.capacity * 100, "%)")
    println(io, "Unique (M): ", t.unique / 1e5)
    println(io, "Occupancy: ", t.unique / t.capacity * 100, "%")
end

function build_mertable!(fn, d::MerTable)
    count_seq = 0
    count_k = 0
    for record in open(fn) |> GzipDecompressorStream |> FASTQ.Reader
        s = FASTQ.sequence(record)
        for m in each(DNAMer{31}, s)
            c = convert(ktype, fwmer(m))
            @inbounds record!(d, c)
        end
        count_seq += 1
        # if count_seq > (1e6)  ## Limit the number of reads for initial testing
        #     break
        # end
    end
    println("Read ", count_seq, " sequences.")
end
