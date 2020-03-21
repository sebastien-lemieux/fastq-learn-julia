using Distributed

addprocs(8)
nprocs()

data = zeros(Int, 1000000000)

@everywhere function addRandom(d::Array{Int}, n)
    println("Starting")
    for i in 1:n
        index = rand(UInt32) % length(d) + 1
        d[index] += 1
    end
    println("Ending")
end

@time for i in 1:8
    r = @spawn addRandom(data, 100000000)
end

12 + 13
