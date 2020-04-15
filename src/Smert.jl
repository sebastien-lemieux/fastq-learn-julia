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
    n = t.unique
    println("Smert: nb. of unique: ", n)
    s = Smert(Array{ktype}(undef, n), Array{ktype}(undef, n), 0,
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
            s.key[i_s] = t.data[i]
            s.count[i_s] = tmp_c
            i_s += 1
        end
    end
    s.n_kc = i_s - 1
    n = s.n_kc

    println("Number of entries passing threshold: ", Int(s.n_kc))

    # Copy only the keys that were allocated.
    s.key = s.key[1:n]
    s.count = s.count[1:n]

    sort(s)
    init_index(s)

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
        s.index[tk] = SmertRange(min(s.index[tk].low, i), max(s.index[tk].high, i))
    end
end

function find_index(t::Smert, k::ktype)
    low = t.index[get_index(k, t)].low
    high = t.index[get_index(k, t)].high
    while high > low
        println((low, high))
        mid = div(low + high, 2)
        if k <= t.key[mid]
            high = mid
        else
            low = mid + 1
        end
    end
    return low
end

function Smert(a::Smert, b::Smert)
    n_a = length(a.key)
    n_b = length(b.key)

    # Empty and to large (will shrink later)
    s = Smert(Array{ktype}(undef, n_a + n_b), Array{ktype}(undef, n_a, n_b), 0,
              fill(SmertRange(typemax(Int), typemin(Int)), 4^indexwidth),
              (31 - indexwidth) * 2)

    i_a = i_b = i_s = 1
    pick(x, i_x) = (s.key[i_s] = x.key[i_x] ; s.count[i_s] = x.count[i_x] ; i_s )
    while i_a <= n_a && i_b <= n_b
        if i_a > n_a
            #pick b
            s.key[i_s] = b.key[i_b]
            s.count[i_s] = b.count[i_b]
            i_s += 1
            i_b += 1
        elseif i_b > n_b
            #pick a
        elseif a.key[i_a] < b.key[i_b]
            #pick a
        elseif b.key[i_b] < a.key[i_a]
            #pick b
        else # they are equal
            #merge a,b
        end
    end

    # k_n = length(t.key)
    # u_n = length(unique_of)
    # nk_n = k_n + u_n
    # new_key = Array{ktype}(undef, nk_n)
    # new_count = Array{ctype}(undef, nk_n)
    # new_index = fill(SmertRange(typemax(Int), typemin(Int)), length(t.index))
    #
    # u_i = k_i = nk_i = 1
    # while (k_i <= k_n || u_i <= u_n)
    #     # println((u_i, k_i, nk_i))
    #     if (k_i <= k_n) && ((u_i > u_n) || (t.key[k_i] < unique_of[u_i]))
    #         new_key[nk_i] = t.key[k_i]
    #         new_count[nk_i] = t.count[k_i]
    #         update_range!(new_index, (new_key[nk_i] >> t.indexshift) + 1, nk_i)
    #         nk_i += 1
    #         k_i += 1
    #     elseif (u_i <= u_n) && ((k_i > k_n) || (unique_of[u_i] < t.key[k_i]))
    #         new_key[nk_i] = unique_of[u_i]
    #         new_count[nk_i] = 0 # Will be added later
    #         update_range!(new_index, (new_key[nk_i] >> t.indexshift) + 1, nk_i)
    #         nk_i += 1
    #         u_i += 1
    #     end
    # end
    # t.key = new_key
    # t.count = new_count
    # t.n_kc = length(new_key)
    # t.index = new_index
    #
    # for i in 1:t.n_of
    #     record!(t, t.overflow[i])
    # end
    # t.n_of = 0
end
