
## Aim for 28.8 sec for a 1M reads. This will give 15m for 250M reads on 8 threads.

include("src/MerTable.jl")

cd("C:\\Users\\slemi\\prog\\Load")

fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"

d = MerTable{UInt64}(Int(1e8))
@time build_mertable!(fn, d)
print(d)
