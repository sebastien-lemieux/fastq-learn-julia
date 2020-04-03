# k-mer table in julia

## Data and objective

__Data__: Raw FASTQ files from RNA-Seq of AML (Acute Myeloid Leukemia) samples. Eventually some
sample-specific attributes to be associated with raw reads.

__Objective__: Demonstrate k-mer level, neural network based prediction of a sample attribute.

__Constraints__:
1. Avoid sample level normalization by using, as primary data source, comparisons between k-mer abundance. Or any value that would not be affected by sequencing depth.
** Rationale: Since it is difficult to demonstrate that a normalization scheme is accurate, it would be important to first demonstrate that it is **needed**.
1. Avoid producing an intermediate "digest" for a sample, prediction needs to proceed from the k-mer table.
** Rationale: This is to limit the amount of computations involved in transforming the k-mer table. Without this constraint, the "digest" could precomputed labels to be predicted.
1. Be robust to missing data. Either missing k-mer sequence (e.g. due to a patient-specific variant) or missing counts (e.g. due to low sequencing depth).
1. Take advantage of BioJulia and Flux

## Roadmap

### Data handling
- [X] Build a k-mer table from a set of FASTQ files for a sample.
- [X] (Opt) Optimize previous step to fit in 10-15 minutes for 250M reads on 8 cores.
- [ ] Manage to fit 500-100 final k-mer tables in 8-16 GB.
- [ ] Load and save the final cohort-wide k-mer matrix.
- [ ] Identify a first sample attribute to predict (probably identification of a relatively easy sub-type... NPM1 mutation?). Obtain, encode, load, save.

### Personal training!
* Learn basic Flux.
* Configure Flux for GPU computing (Windows [early dev] + Linux [Later res. gathering])

### Main task
* Research an algo to get a question-driven prediction (request a variable number of k-mer count comparisons). Look into reinforcement learning strategy.
** Alternatively, identify a small subset of the comparisons that would be sufficient to drive prediction (a la L1000). Start with random selection, then random + prioritize.

### Reporting
- [X] Preparation of a time sheet to track time investment.
- [X] Setup basic Github repository structure and establish initial workflow from Atom/Juno.
- Final report (where, structure, res. production, fig. prep...)

Time sheet: https://docs.google.com/spreadsheets/d/1LEfk4cWtw-cFPHTBzwMCpgf5g2ybY8MZ8gd2DBbS5bo/edit?usp=sharing
