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

function count_total(t::Smert)
    sum = 0
    for i in 1:t.n_kc
        sum += t.count[i]
    end
    return sum
end


function Smert(a::Smert, b::Smert)
    n_a = length(a.key)
    n_b = length(b.key)

    # Empty and too large (will shrink later)
    s = Smert(Array{ktype}(undef, n_a + n_b), Array{ktype}(undef, n_a + n_b), 0,
              fill(SmertRange(typemax(Int), typemin(Int)), length(a.index)),
              a.indexshift)

    i_a = i_b = i_s = 1

    while i_a <= n_a && i_b <= n_b
        # if i_s % 1000 == 0
        #     println(i_s)
        # end
        if i_a > n_a
            #pick b
            s.key[i_s] = b.key[i_b]
            s.count[i_s] = b.count[i_b]
            i_s += 1
            i_b += 1
        elseif i_b > n_b
            #pick a
            s.key[i_s] = a.key[i_a]
            s.count[i_s] = a.count[i_a]
            i_s += 1
            i_a += 1
        elseif a.key[i_a] < b.key[i_b]
            #pick a
            s.key[i_s] = a.key[i_a]
            s.count[i_s] = a.count[i_a]
            i_s += 1
            i_a += 1
        elseif b.key[i_b] < a.key[i_a]
            #pick b
            s.key[i_s] = b.key[i_b]
            s.count[i_s] = b.count[i_b]
            i_s += 1
            i_b += 1
        else # they are equal
            #merge a,b
            s.key[i_s] = a.key[i_a]
            s.count[i_s] = a.count[i_a] + b.count[i_b]
            i_s += 1
            i_a += 1
            i_b += 1
        end
    end

    s.n_kc = i_s - 1
    s.key = s.key[1:(s.n_kc)]
    s.count = s.count[1:(s.n_kc)]

    println(s.n_kc)

    # s is already sorted, only needs indexing
    init_index(s)

    return s
end
