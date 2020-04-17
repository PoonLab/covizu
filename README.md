## Objectives

* To develop an open-source pipeline for the alignment and phylogenetic analysis of SARS-COV-2 genome sequences from the GISAID database.
* To experiment with new visualization methods on the outputs of these analyses.
* To attempt to extract actionable information for public health responses to the ongoing pandemic.

## Current workflow
1. Sequences are bulk downloaded from the GISAID database.  All developers have signed the GISAID data access agreement, and sequences are not being re-distributed.
2. Sequences are aligned pairwise against the SARS-COV-2 reference genome using the Procrustean method implemented in [gotoh2](http://github.com/ArtPoon/gotoh2).  This module provides a method that progressively updates an existing alignment file with new sequence records, avoiding the re-alignment of previously released genomes.
3. Sequences are filtered using `filtering.py` for entries that are derived from non-human sources, incomplete genomes, and genomes that contain >5% fully ambiguous base calls (`N`s).
4. A pairwise genetic distance matrix is generated using [TN93](http://github.com/veg/tn93).
5. Clusters of effectively identical genome sequences are identified using the Python script `clustering.py` and a minimum spanning tree is generated using the [networkx](https://networkx.github.io/) module.  The tree is partitioned into clusters by cutting branches at nodes with a high degree size.
6. A plot is generated from a given cluster using the R script `draw-mst.R`.
