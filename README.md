# h3n2-segment-rates

## Segment-specific evolutionary rates in human seasonal H3N2 influenza A: linking molecular clock heterogeneity to positive selection on surface glycoproteins.

## Question

Do the eight influenza A genome segments evolve at different rates in human seasonal H3N2, and does the surface-protein vs. internal-protein distinction explain the difference?

Hypothesis. HA (segment 4) and NA (segment 6) encode the surface glycoproteins exposed to host antibodies. They should show both higher substitution rates and more positively selected sites than the six internal segments (PB2, PB1, PA, NP, M, NS).

Why it matters. Rate heterogeneity across segments is the mechanism behind annual vaccine reformulation: the segments the immune system sees are the ones that change fastest, which is why vaccine strain prediction is hard.

## Study design

Taxa ~150 H3N2 strains (target 156 = 6/year × 26 years)
Loci 8 — one per genome segment
Time span 2000–2025, evenly subsampled
Host Homo sapiens only
Completeness All 8 segments required per strain
Database NCBI Virus (GenBank)

## Getting the raw data

Source: NCBI Virus — https://www.ncbi.nlm.nih.gov/labs/virus/vssi/

## Repository structure

h3n2-segment-rates/
├── README.md
├── ROADMAP.md
├── data/
│ ├── raw/ # NCBI download (gitignored — large)
│ └── processed/ # strain-matched, subsampled FASTA + metadata
├── scripts/
│ ├── 00_validate_download.R
│ ├── 01_match_subsample.R
│ ├── 02_align.sh
│ ├── 03_modelfinder_trees.sh
│ ├── 04_mrbayes.sh
......
├── results/  
├── figures/
├── paper.Rmd
└── references.bib

## Software

Tool Purpose
MAFFT per-segment alignment
IQ-TREE 2 ModelFinder, ML trees, UFBoot + SH-aLRT, concordance factors
MrBayes Bayesian tree (reduced taxon set, method comparison)
TreeTime tip-dated molecular clock, per-segment rates
HyPhy BUSTED, MEME, FUBAR selection tests
