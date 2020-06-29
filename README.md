# Genomic Surveillance Reveals Multiple Introductions of SARS-CoV-2 into Northern California 

**Xianding Deng**, **Wei Gu**, **Scot Federman**, **Louis du Plessis**, Oliver G. Pybus, Nuno Faria, Candace Wang, Guixia Yu, Chao-Yang Pan, Hugo Guevara, Alicia Sotomayor-Gonzalez, Kelsey Zorn, Allan Gopez, Venice Servellita, Elaine Hsu, Steve Miller, Trevor Bedford, Alexander L. Greninger, Pavitra Roychoudhury, Lea M. Starita, Michael Famulare, Helen Y. Chu, Jay Shendure, Keith R. Jerome, Catie Anderson, Karthik Gangavarapu, Mark Zeller, Emily Spencer, Kristian G. Andersen, Duncan MacCannell, Clinton R. Paden, Yan Li, Jing Zhang, Suxiang Tong, Gregory Armstrong, Scott Morrow, Matthew Willis, Bela T. Matyas, Sundari Mase, Olivia Kasirye, Maggie Park, Curtis Chan, Alexander T. Yu, Shua J. Chai, Elsa Villarino, Brandon Bonin, Debra A. Wadford, and **Charles Y. Chiu**

[![DOI](https://zenodo.org/badge/259707507.svg)](https://zenodo.org/badge/latestdoi/259707507)

---

This repository contains the data files, scripts and workflows necessary to reproduce the phylogenetic analyses and figures presented in **Deng et al., Science, 2020** ([https://doi.org/10.1126/science.abb9263](https://doi.org/10.1126/science.abb9263)). Some of the scripts may need some adjustment depending on the local setup. 

Note that because of the GISAID [terms of use](https://www.gisaid.org/registration/terms-of-use/) genomic sequences cannot be shared in this repository. Instead, we make the GISAID accessions available.


## Abstract

_The novel coronavirus SARS-CoV-2 has caused a major global pandemic. Here we investigate the genomic epidemiology of SARS-CoV-2 in Northern California using samples from returning travelers, cruise ship passengers, and cases of community transmission. Virus genomes were recovered from 36 SARS-CoV-2 positive patients from the end of January through March 20th. Phylogenetic analyses revealed at least 7 different SARS-CoV-2 lineages, suggesting multiple independent introductions of the virus into the state. Virus genomes from passengers of the Grand Princess cruise ship clustered with those from an established epidemic in Washington State, including the first reported case in the United States. We also detected evidence for presumptive transmission of different lineages between communities. These findings support the universal implementation of widespread testing, contact tracing, social distancing, and travel restrictions to mitigate the virus spread._


## Dependencies

- QuickTree
- PhyML v3.3
- BioPython
- ggplot2, ggtree, ape, phytools, treeio, lubridate, cowplot
- [beastio](https://github.com/laduplessis/beastio): Commit #b18caa6


## Data

- [`gisaid_cov2020_CA_789seqs.MSA.trimmed_taxa.csv`](https://github.com/laduplessis/SARS-CoV-2_Northern_California_genomic_surveillance/blob/master/data/gisaid_cov2020_CA_789seqs.MSA.trimmed_taxa.csv): Taxon labels of sequences in the global alignment (with trimmed ends).
- [`WA1_taxa.csv`](https://github.com/laduplessis/SARS-CoV-2_Northern_California_genomic_surveillance/blob/master/data/WA1_taxa.csv): Taxon labels of sequences in the WA1 cluster alignment 
- [`MN908947.3.fasta`](https://github.com/laduplessis/SARS-CoV-2_Northern_California_genomic_surveillance/blob/master/data/MN908947.3.fasta): This is the reference sequence (Wuhan-Hu-1)
- [`california_clusters.csv`](https://github.com/laduplessis/SARS-CoV-2_Northern_California_genomic_surveillance/blob/master/data/california_clusters.csv): Sample metadata and definitions of the clusters.


## Phylogenetic analyses

### Prepare alignments 

```bash
# Rename ids so it can be more easily parsed in files
tr '/' '_' < data/gisaid_cov2020_CA_789seqs.MSA.trimmed.fasta > results/alignments/gisaid_cov2020_CA_789seqs.MSA.trimmed.fasta

# Create USA only alignment
grep USA -A 1 results/alignments/gisaid_cov2020_CA_789seqs.MSA.trimmed.fasta > results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.fasta

# Remove ambiguous sites 
python scripts/convertmsa.py -i results/alignments/gisaid_cov2020_CA_789seqs.MSA.trimmed.fasta -o results/alignments/gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.fasta -u -r
python scripts/convertmsa.py -i results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.fasta -o results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.noAmbiguous.fasta -u -r

```
Save alignments as Phylip files with full and padded names using Aliview.


### Sequence statistics for alignment figures

Remove all gaps and indels in the USA only alignment so it has the same length as the reference sequence (29,903bp)

```bash
python scripts/msastats.py -i results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.fasta -r data/MN908947.3.fasta -f "|" -s 0 -o results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed/
python scripts/msastats.py -i results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.noAmbiguous.fasta -r data/MN908947.3.fasta -f "|" -s 0 -o results/alignments/gisaid_cov2020_CA_152seqs.MSA.trimmed.noAmbiguous/

```

### PhyML (SH-like branch supports)

Approximate likelihood ratio test returning SH-like non-parametric branch supports.

```bash
# HKY
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -a 1 --quiet --run_id HKY > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy --quiet --run_id HKY+G > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+G.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -v e --quiet --run_id HKY+G+I > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+G+I.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -a 1 -v e --quiet --run_id HKY+I > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+I.out &

# HKY (no ambiguous sites)
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -a 1 --quiet --run_id HKY > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy --quiet --run_id HKY+G > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+G.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -v e --quiet --run_id HKY+G+I > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+G+I.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -a 1 -v e --quiet --run_id HKY+I > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+I.out &

```


### PhyML (aLRT Chi2)

Approximate likelihood ratio test returning Chi2-based parametric branch supports.

```bash
# HKY
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -b -2 -a 1 --quiet --run_id HKY > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -b -2 --quiet --run_id HKY+G > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+G.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -b -2 -v e --quiet --run_id HKY+G+I > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+G+I.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.phy -b -2 -a 1 -v e --quiet --run_id HKY+I > gisaid_cov2020_CA_789seqs.MSA.trimmed_HKY+I.out &


# HKY (no ambiguous sites)
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -b -2 -a 1 --quiet --run_id HKY > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -b -2 --quiet --run_id HKY+G > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+G.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -b -2 -v e --quiet --run_id HKY+G+I > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+G+I.out &
~/phyml -i gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous.phy -b -2 -a 1 -v e --quiet --run_id HKY+I > gisaid_cov2020_CA_789seqs.MSA.trimmed.noAmbiguous_HKY+I.out &

```

### WA1 only alignment

Create alignments

```bash
mkdir results/WA1/
touch results/WA1/WA1.fasta
for SEQ in `cat data/2020-04-27/WA1_taxa.csv`
do 
	grep -A 1 $SEQ results/alignments/gisaid_cov2020_CA_789seqs.MSA.trimmed.fasta >> results/WA1/WA1.fasta
done;

# Remove ambiguous sites and get base composition of each column in the alignment
python scripts/convertmsa.py -i results/WA1/WA1.fasta -o results/WA1/WA1.noAmbiguous.fasta -u -r
python scripts/msastats.py -i results/WA1/WA1.noAmbiguous.fasta -f "|" -s 0 -o results/WA1/WA1_stats/

# Maximum of 4 sequences sharing an ambiguous site/gap (delete 528 sites)
python scripts/trimsequences.py -a results/WA1/WA1.noAmbiguous.fasta -H results/WA1/WA1_stats/WA1.noAmbiguous.hist.csv -o results/WA1/ -p WA1.noAmbiguous.0.05 -c 0.05

# No shared ambiguous sites/gaps (delete 4773 sites)
python scripts/trimsequences.py -a results/WA1/WA1.noAmbiguous.fasta -H results/WA1/WA1_stats/WA1.noAmbiguous.hist.csv -o results/WA1/ -p WA1.noAmbiguous.0.02 -c 0.02

# No ambiguous sites/gaps at all (delete 13521 sites)
python scripts/trimsequences.py -a results/WA1/WA1.noAmbiguous.fasta -H results/WA1/WA1_stats/WA1.noAmbiguous.hist.csv -o results/WA1/ -p WA1.noAmbiguous.0 -c 0

```

Save the alignments as Phylip files and run below for the maximum-likelihood trees

```bash
phyml -i WA1.noAmbiguous.0.phy -b -2 --quiet --run_id HKY+G > WA1.noAmbiguous.0.HKY+G.out &
phyml -i WA1.noAmbiguous.0.02.phy -b -2 --quiet --run_id HKY+G > WA1.noAmbiguous.0.02.HKY+G.out &
phyml -i WA1.noAmbiguous.0.05.phy -b -2 --quiet --run_id HKY+G > WA1.noAmbiguous.0.05.HKY+G.out &

```

## Figures

To reproduce the figures, run the notebook `scripts/Figures.Rmd`, which creates [Figures.pdf](https://github.com/laduplessis/SARS-CoV-2_/blob/master/scripts/Figures.pdf). Note that the notebook requires the SNP and gap statistics of the alignment for some figures, which are not stored on the repository (because of GISAID terms of use) and that the figures require some post-processing.
