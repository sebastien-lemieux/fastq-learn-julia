
## Aim for 28.8 sec for a 1M reads. This will give 15m for 250M reads on 8 threads.

include("src/MerTable.jl")
include("src/Smert.jl")

cd("C:\\Users\\slemi\\prog\\Load")

# function build_mertable(directory::String, d::MerTable{E}) where {E}
#     all_files = readdir(directory)
#     Threads.@threads for f in all_files
#         println("Starting ", f)
#         build_mertable!(directory * f, d)
#         println("Ending   ", f)
#     end
#     return d
# end


fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"
# d = Smert(10000000, 13)
# @time build_smert!(fn, d)

t = MerTable(1e8)
@time build_mertable!(fn, t)  # 15-18s/M
@time s = Smert(t, 12, 0) # 24s/M (fixed time?)
print(t)
