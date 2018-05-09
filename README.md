# SNV filters
Somatypus post-processing script (in R), after removing variants with flags 'badReads', 'MQ', 'strandBias', 'SC', 'QD'.

<b>1. Homozygous-in-reference filter.</b>
SNVs called with VAF > 0.9 in reference sample are discarded.

<b>2. Strand bias filter.</b>
In regions with total coverage ≥ 11 across all samples, we reject variant calls with less than 20 % support on either forward or reverse sequencing strands. In regions with total coverage < 11 reads across all samples, we remove variants that had less than two supporting reads in either forward or reverse direction.

<b>3. Simple repeats regions filter.</b>
SNVs lying within a 5 BP window around a simple repeat regions, as annotated by Tandem Repeat Finder, are discarded.

<b>4. Regions filter.</b>
The Tasmanian devil reference genome (Devil7.1) is a scaffold-level assembly, consisting of 237,291 contigs assembled into 35,974 scaffolds. We reject any variant mapping within 500 BP from the start or end of a contig, or within 1000 BP from the start of end of a scaffold.
