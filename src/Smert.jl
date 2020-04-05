using FASTX
using CodecZlib
using BioSequences

include("common.jl")

struct SmertRange
    low::Int
    high::Int
end

mutable struct Smert # Sorted Mer Table
    key::Array{ktype}
    count::Array{ctype}
    n_kc::ctype  # nb. of keys in key-count table
    index::Array{SmertRange}
    indexshift::ktype
end

function Smert(t::MerTable, indexwidth, subtcount = 0)
    println("yeah")
    n = t.unique
    println("Smert: nb. of unique: ", n)
    # s = Smert(Array{ktype}(undef, n), Array{ktype}(undef, n), 0,
    #           fill(SmertRange(typemax(Int), typemin(Int)), 4^indexwidth),
    #           (31 - indexwidth) * 2)

    s = Smert(zeros(ktype, n), zeros(ctype, n), 0,
              fill(SmertRange(typemax(Int), typemin(Int)), 4^indexwidth),
              (31 - indexwidth) * 2)
    println("index size: ", 4^indexwidth)

    i_s = 1
    for i in 1:length(t.data)
        tmp_c = t.count[i]
        if tmp_c > subtcount
            tmp_c -= subtcount
        else
            tmp_c = 0
        end
        if tmp_c > 0
            if i_s > n
                println("yurk!!!!!")
                println("yurk!!!!!")
                println("yurk!!!!!")
                println("yurk!!!!!")
                println("yurk!!!!!")
                break
            end
            # s.key[i_s] = t.data[i]
            # s.count[i_s] = tmp_c
            i_s += 1
        end
    end
    s.n_kc = i_s - 1
    # n = s.n_kc

    println("Number of entries passing threshold: ", Int(s.n_kc))

    # Copy only the keys that were allocated.
    # s.key = s.key[1:n]
    # s.count = s.count[1:n]

    # sort(s)
    # create_index(s)

    return s
end

import Base.sort  ## needed to extend an existing function
function sort(s::Smert)
    println("Sorting")
    p = sortperm(s.key)
    s.key = s.key[p]
    s.count = s.count[p]
end

get_index(key::ktype, t::Smert) = key >> t.indexshift + 1

function init_index(s::Smert)
    println("Indexing")
    for i in 1:length(s.key)
        tk = get_index(s.key[i], s)
        s.index[tk].low = min(s.index[tk].low, i)
        s.index[tk].high = min(s.index[tk].high, i)
    end
end

function find_index(t::Smert, k::ktype)
    @inbounds low = t.index[get_index(k, t)].low
    @inbounds high = t.index[get_index(k, t)].high
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

# function update_range!(index::Array{SmertRange}, tk::ktype, nk_i)
#     @inbounds index[tk] = SmertRange(min(index[tk].low, nk_i),
#                            max(index[tk].high, nk_i))
# end

# function merge_of(t::Smert)
#     unique_of = unique(sort(t.overflow[1:t.n_of]))
#     k_n = length(t.key)
#     u_n = length(unique_of)
#     nk_n = k_n + u_n
#     new_key = Array{ktype}(undef, nk_n)
#     new_count = Array{ctype}(undef, nk_n)
#     new_index = fill(SmertRange(typemax(Int), typemin(Int)), length(t.index))
#
#     u_i = k_i = nk_i = 1
#     while (k_i <= k_n || u_i <= u_n)
#         # println((u_i, k_i, nk_i))
#         if (k_i <= k_n) && ((u_i > u_n) || (t.key[k_i] < unique_of[u_i]))
#             new_key[nk_i] = t.key[k_i]
#             new_count[nk_i] = t.count[k_i]
#             update_range!(new_index, (new_key[nk_i] >> t.indexshift) + 1, nk_i)
#             nk_i += 1
#             k_i += 1
#         elseif (u_i <= u_n) && ((k_i > k_n) || (unique_of[u_i] < t.key[k_i]))
#             new_key[nk_i] = unique_of[u_i]
#             new_count[nk_i] = 0 # Will be added later
#             update_range!(new_index, (new_key[nk_i] >> t.indexshift) + 1, nk_i)
#             nk_i += 1
#             u_i += 1
#         end
#     end
#     t.key = new_key
#     t.count = new_count
#     t.n_kc = length(new_key)
#     t.index = new_index
#
#     for i in 1:t.n_of
#         record!(t, t.overflow[i])
#     end
#     t.n_of = 0
# end

# function record!(t::Smert, m::ktype)
#     if t.n_kc > 0
#         i = find_index(t, m)
#         @inbounds if i != typemax(Int) && t.key[i] == m
#             t.count[i] += 1
#             return (t.n_kc, t.n_of)
#         end
#     end
#     t.n_of += 1
#     t.overflow[t.n_of] = m
#     if t.n_of == length(t.overflow)
#         println("Auto merge")
#         merge_of(t)
#     end
#     return (t.n_kc, t.n_of)
# end
