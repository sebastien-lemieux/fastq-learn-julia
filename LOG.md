2020-03-05 (J01, 3h):
* 1800-2100
* Install on windows (Juno, Julia, Atom, BioJulia, FASTX, CodecZlib)
* Load fastq files (FASTX, BioJulia)
* Read on biojulia benchmark comparison with Seq language
* Load compressed fastq files
* Multi-thread version of the loader (8-files in parallel)
* cpu-bound task
* overclock to 4.8GHz (temp stable, fans low)
* On 8 threads, reaches 3.6M reads per sec (est. 1m for a whole dataset, 224M reads)
* Next: encode in biosequences, test amount of ram to keep sequences, then extract k-mers

2020-03-06 (J02, 1h45):
* Extract sequence from record, create biosequence
* First successful attempt at building k-mer table (4M, k=5, 81s)
* Most time spent hashing
* Using convert as a specialized hash: (13.1s)
* Appears much longer with k=14 (never finishes???), maybe try to rebuild from scratch = understand + language

2020-03-09 (J05, 1h30):
* Windows update (yuk!) + Read on stream
* Basic struct to create an iterable (iter_test.jl)
* Basic fastq parser as an iterable (load3.jl)

2020-03-12 (J08, 1h10):
* Read on parametric types, not entirely clear... Not using it.
* More complex iterable to extract k-mer, currently as string
* Next: convert MerGenerator to use a UInt64 mer representation, need DNA convert to 2 bits.

2020-03-15 (J11, 2h15):
* Finished a "from scratch" 2-bit implementation of k-mer table (N -> A currently).
* Performance down to 18k reads per sec! (for any k up to 31).
* Re-test of the biojulia implementation, time increase dramatically around k = 11. Due to the dict.
* With k=31 the biojulia impl. without storing in the Dict reaches 423k reads per sec.
* Next: re-read (https://discourse.julialang.org/t/developing-a-kmer-counter-in-julia/292/2)+
  implement a basic hash table for biojulia dnamer?

2020-03-16 (J12, 2h):
* Initial github setup, explore the use of github within Atom/juno
* Expanded the definition of MerTable.
* Experimentation with parametric types (now working!).
* Primitive hash table (with collision).
* 6.8 s/M (without checks).  13.7 with checks.
* Next: Resolve collision.

2020-03-18 (J14, 2h):
* Added linear probing to resolve collisions. 12.0 s/M.
* Explored optimization: no gains from @inbounds, most time spent in correct mapping.
* Little expected gains from anything more fancy than linear probing.
* Not ideal reporting of collision frequency.
* Naive parallel implementation (incorrect), issue with GC / Mem allocation.
  Need better understanding of parallel (threads) vs. mutable arrays.
* Next: check @spawn + Task

2020-03-21 (J17, 4h):
* Reading on multiprocessing using Distributed module... need to test mem. management.
  Confirmed that arrays are copied when @spawn calling.
* Started work on an alternative MerTable (branch sorted_table).
* New approach is very memory efficient, but very slow (~105 s/M!).
* Add 12nt index to the sorted structure (a la STAR).

2020-03-24 (2h):
* First attempt at a 12-bp preindex, 50.8 s/M.
* Multiple attempt (failed) at profiling, including allocation profiling.
* Allocations can be tracked using option: --track-allocation=all (juno:settings)
* Found out that most allocation are done when calling getindex() on an Array.
  Difficult to avoid. This fits with most lines taking time in the profiling
  being [] accesses on arrays. @inbounds seems to reduce some allocations.

2020-03-27 (1h30):
* Need to combine both approach. Use hash to build the table then transfer to a
  sorted table to save space to hold the whole cohort.
* Lots need to be done to complete the merge: sort the Smert, create index, query
  the MerTable, handle removal of low count Mers...

2020-04-03 (1h40):
* Update to Readme mostly. Some performance testing.

2020-04-05 (1h45):
* Update installation to Julia 1.4.0
* Started work on the MerTable import, stuck on a bug (missing string in println!) in Smert
  constructor.

2020-04-15 (11:45-14:15, 15h- ):
* Download data for a new sample (complete): 18H175
* Intro to Windows' powershell, CLI ssh/scp, install of gzip
* The strange bug (04-05) is actually in Atom/juno function to send a selected lines to the repl.
* Smert:index built (untested).
* Now loads two files in two Smert, ready to work on merge.
