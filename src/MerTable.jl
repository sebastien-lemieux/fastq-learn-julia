using FASTX
using CodecZlib
using BioSequences

ktype = UInt64    # Up to 32-mer
ctype = UInt32  # Up to 4G occurrences

struct TableRange
    low::Int
    high::Int
end

mutable struct MerTable
    key::Array{ktype}
    count::Array{ctype}
    overflow::Array{ktype}
    n_kc::ctype  # nb. of keys in key-count table
    n_of::ctype  # nb. of keys in overflow
    index::Array{TableRange}
end

function MerTable(overflowsize)
    return MerTable(Array{ktype}(undef, 0), Array{ktype}(undef, 0),
                    Array{ktype}(undef, UInt(overflowsize)),
                    0, 0,
                    Array{TableRange}(undef, 0))
end

function find_index(t::MerTable, k::ktype)
    low = 1
    high = Int(length(t.key))
    while high > low
        # println((low, high))
        mid = div(low + high, 2)
        if k <= t.key[mid]
            high = mid
        else
            low = mid + 1
        end
    end
    return low
end

# d = MerTable(10)
# d.key = [1, 2, 3, 4, 5, 6, 7]
# find_index(d, ktype(6))

# function find_index(t::MerTable, k::ktype, low::Int=1, high::Int=length(t.key))
#     # println((low, high))
#     if low == high
#         return low
#     end
#     mid = div(low + high, 2)
#     if k <= t.key[mid]
#         return find_index(t, k, low, mid)
#     else
#         return find_index(t, k, mid + 1, high)
#     end
# end

function merge_of(t::MerTable)
    unique_of = unique(sort(t.overflow[1:t.n_of]))
    k_n = length(t.key)
    u_n = length(unique_of)
    nk_n = k_n + u_n
    new_key = Array{ktype}(undef, nk_n)
    new_count = Array{ctype}(undef, nk_n)

    **** prepare the index ****

    u_i = k_i = nk_i = 1
    while (k_i <= k_n || u_i <= u_n)
        # println((u_i, k_i, nk_i))
        if (k_i <= k_n) && ((u_i > u_n) || (t.key[k_i] < unique_of[u_i]))
            new_key[nk_i] = t.key[k_i]
            new_count[nk_i] = t.count[k_i]
            nk_i += 1
            k_i += 1
        elseif (u_i <= u_n) && ((k_i > k_n) || (unique_of[u_i] < t.key[k_i]))
            new_key[nk_i] = unique_of[u_i]
            new_count[nk_i] = 0 # Will be added later
            nk_i += 1
            u_i += 1
        end
    end
    t.key = new_key
    t.count = new_count
    t.n_kc = length(new_key)

    for i in 1:t.n_of
        record!(t, t.overflow[i])
    end
    t.n_of = 0
end

function record!(t::MerTable, m::ktype)
    if t.n_kc > 0
        i = find_index(t, m)
        if t.key[i] == m
            t.count[i] += 1
            return (t.n_kc, t.n_of)
        end
    end
    t.n_of += 1
    t.overflow[t.n_of] = m
    if t.n_of == length(t.overflow)
        println("Auto merge")
        merge_of(t)
    end
    return (t.n_kc, t.n_of)
end

# m = MerTable(10000)
# for i in (1:25)
#     record!(m, ktype(i))
# end
# merge_of(m)
# for i in (1:12)
#     println(ktype(i + 30))
#     record!(m, ktype(i + 30))
# end

# @time for i in 1:1000000
#     record!(m, ktype(rand(1:10000000)))
# end


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
        if count_seq % 1000 == 0
            println(count_seq)
        end
        if count_seq > (2e5)
            break
        end
    end
end
