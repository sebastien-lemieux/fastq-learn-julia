using GZip

function nuctonum(nt)
    if nt == 'A'
        return 0
    elseif nt == 'C'
        return 1
    elseif nt == 'G'
        return 2
    elseif nt == 'T'
        return 3
    else
        return 4
    end
end

function treat_seq!(d, k, seq)
    curr = 0
    curr_n = 0
    mask = 1 << (2 * k) - 1
    for i = 1:length(seq)
        num = nuctonum(seq[i])
        if num == 4
            curr = 0
            curr_n = 0
        else
            curr = ((curr << 2) + num) & mask
            curr_n += 1
            if curr_n >= k
                d[curr+1] += 1
            end
        end
    end
end

cd("/Users/slemi/prog")

@time GZip.open("10H005/10H005_ACTTGA_L002_R1_004.fastq.gz", "r") do f

    k = 14
    d = zeros(Int32, 4^k)

    count = 0
    readline(f)       # title
    t = @elapsed for seq in eachline(f)
        readline(f)       # +
        readline(f)       # quality
        readline(f)       # next title
        count = count + 1
        treat_seq!(d, k, seq)
        # if count == 400000
        #     break
        # end
    end
    #print(d)
    print("$count lines read in $t seconds ($(count / t) l/s)")
end
