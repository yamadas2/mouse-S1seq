# mouse S1-seq
This pipeline contains code used for read processing and mapping of mouse S1-seq data. Example R scripts to visualize S1-seq profiles are also included. The pipeline can be used to map reads of S1-seq and variants such as ExoT-seq and Exo7/T-seq, against references other than the mouse genome.

S1-seq and its variant methods have been developed to visualize the global distribution of endpoints of meiotic DNA double-strand break processing. In those methods, nuclease S1 or other single-strand specific nucleases are used to trim ssDNA tails at resected DSB ends, which allows ligation of biotinylated sequencing adapters to resected DSB ends and mapping the ends of the resection trancts.

The pipeline starts by trimming paired-end sequence reads, followed by mapping onto a reference genome. After removal of duplicated reads, it counts uniquely and properly mapped reads at the nucleotide next to the biotinylated adapter was mapped (this corresponds to the resection endpoint). 

Example R scripts included are for generating mean profiles and heatmaps around mouse SPO11-oligo hotspot centers, calculating mean resection lengths and scatter plots to check correlation with SPO11-oligo maps and reproducibility between maps.

## Installation
Installation instructions are described in setup.txt in mapping/scripts/MouseResection/main.

## Usage
Usage is described in usage.md in mapping/scripts/MouseResection/main.
