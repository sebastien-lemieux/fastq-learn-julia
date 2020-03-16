# n = 8
# t = Array{Int}(undef, n)
# Threads.@threads for i in 1:n
#     t[i] = count_fastq("../10H005/10H005_ACTTGA_L001_R1_$(@sprintf("%03d",i)).fastq.gz")
# end
# println("vlan")
#
# for i in 1:n
#     println(t[i])
# end
