include("src/MerTable.jl")

cd("C:\\Users\\slemi\\prog\\Load")

fn = "../10H005/10H005_ACTTGA_L001_R1_001.fastq.gz"

@time d = build_ktable(fn)
print(d)
