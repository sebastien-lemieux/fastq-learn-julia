2020-03-05 (J001, 3h):
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

2020-03-06 (J002, 9:45-11:30, 1h45, 4h45):
* Extract sequence from record, create biosequence
* First successful attempt at building k-mer table (4M, k=5, 81s)
* Most time spent hashing
* Using convert as a specialized hash: (13.1s)
* Appears much longer with k=14 (never finishes???), maybe try to rebuild from scratch = understand + language

2020-03-09 (J005, 10:15-11:45, 1h30, 6h15):
* Windows update (yuk!) + Read on stream
* Basic struct to create an iterable (iter_test.jl)
* Basic fastq parser as an iterable (load3.jl)

2020-03-12 (J008, 19:50-21:00, 1h10, 7h25):
* Read on parametric types, not entirely clear... Not using it.
* More complex iterable to extract k-mer, currently as string
* Next: convert MerGenerator to use a UInt64 mer representation, need DNA convert to 2 bits.

2020-03-15 (J011, 14:45-17:00, 2h15, 9h40):
* Finished a "from scratch" 2-bit implementation of k-mer table (N -> A currently).
* Performance down to 18k reads per sec! (for any k up to 31).
* Re-test of the biojulia implementation, time increase dramatically around k = 11. Due to the dict.
* With k=31 the biojulia impl. without storing in the Dict reaches 423k reads per sec.
* Next: re-read (https://discourse.julialang.org/t/developing-a-kmer-counter-in-julia/292/2)+
  implement a basic hash table for biojulia dnamer?

2020-03-16 (J012, 15:15-15:40, 17h35-17h40, 17h55-):
* Initial github setup, explore the use of github within Atom/juno
* Expanded the definition of MerTable.
* Experimentation with parametric types (now working!).
* Primitive hash table (with collision).
* 6.8 s/M (without checks).
* Next: Resolve collision.
