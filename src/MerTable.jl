using FASTX
using CodecZlib
using Printf
using BioSequences

# import Base.hash
# function hash(x::DNAMer{31})::UInt64
#     return convert(UInt64,x)
# end

function build_ktable(fn)
    d = Dict{DNAMer{31},Int}()
    count_seq = 0
    count_k = 0
    for record in open(fn) |> GzipDecompressorStream |> FASTQ.Reader
        s = FASTQ.sequence(record)
        for m in each(DNAMer{31}, s)
            c = fwmer(m)
            count_k += haskey(d, c)
            # if haskey(d, c)
            #     d[c] += 1
            # else
            #     d[c] = 1
            # end
        end
        count_seq += 1
        if count_seq > 1e6
            break
        end
    end
    println(count_k)
    return d
end
