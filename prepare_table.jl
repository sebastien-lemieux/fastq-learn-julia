
## Aim for 28.8 sec for a 1M reads. This will give 15m for 250M reads on 8 threads.

include("src/MerTable.jl")
include("src/Smert.jl")

cd("/Users/slemi/prog/data")

# function build_mertable(directory::String, d::MerTable{E}) where {E}
#     all_files = readdir(directory)
#     Threads.@threads for f in all_files
#         println("Starting ", f)
#         build_mertable!(directory * f, d)
#         println("Ending   ", f)
#     end
#     return d
# end


# fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"
fn_1 = "NS.1134.001.NEBNext_dual_i7_E3---NEBNext_dual_i5_E3.18H175_R1.fastq.gz"
fn_2 = "NS.1134.001.NEBNext_dual_i7_E3---NEBNext_dual_i5_E3.18H175_R2.fastq.gz"

t = MerTable(1e9)
@time build_mertable!(fn_1, t)  # 15-18s/M
@time s_1 = Smert(t, 12, 1) # 24s/M (fixed time?)

t = MerTable(1e9)
@time build_mertable!(fn_2, t)  # 15-18s/M
@time s_2 = Smert(t, 12, 1) # 24s/M (fixed time?)



seq = DNAMer{31}("AAAAAAAATTTAATCAAGTGAAACGTAATAA")

find_index(s_1, convert(ktype, seq))
find_index(s_2, convert(ktype, seq))

@time s_m = Smert(s_1, s_2)

find_index(s_m, convert(ktype, seq))
