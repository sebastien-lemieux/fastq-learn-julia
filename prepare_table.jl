
## Aim for 28.8 sec for a 1M reads. This will give 15m for 250M reads on 8 threads.

include("src/MerTable.jl")

cd("C:\\Users\\slemi\\prog\\Load")

fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"

function build_mertable(directory::String, d::MerTable{E}) where {E}
    all_files = readdir(directory)
    Threads.@threads for f in all_files
        println("Starting ", f)
        build_mertable!(directory * f, d)
        println("Ending   ", f)
    end
    return d
end


# d = MerTable{UInt64}(Int(1e8))
# @time build_mertable!(fn, d)

d = 0
GC.gc()
@time d = MerTable{UInt64}(Int(1e9))
@time d = build_mertable("../10H005/", d)


print(d)

# 4 12.5
# 8
