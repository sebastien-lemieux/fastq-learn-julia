using FASTX
using CodecZlib
using Printf
using BioSequences

# function count_fastq(fn)
#     println("starting $(Threads.threadid())")
#     count = 0
#     for record in open(fn) |> GzipDecompressorStream |> FASTQ.Reader
#         count += 1
#     end
#
#     println("ending $(Threads.threadid())")
#     return(count)
# end

# 128s
cd("C:\\Users\\slemi\\prog\\Load")

import Base.hash
function hash(x::DNAMer{31})::UInt64
    return convert(UInt64,x)
end

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

fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"
# @profiler d = build_ktable(fn)
@time d = build_ktable(fn)
print(d)

# 460k reads per sec, 1 thread

# n = 8
# t = Array{Int}(undef, n)
# Threads.@threads for i in 1:n
#     t[i] = count_fastq("../10H005/10H005_ACTTGA_L001_R1_$(@sprintf("%03d",i)).fastq.gz")
# end
# println("vlan")
#
# for i in 1:n
#     println(t[i])
# end
