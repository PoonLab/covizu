## Objectives

* To develop an open-source pipeline for the alignment and phylogenetic analysis of SARS-COV-2 genome sequences from the GISAID database.
* To experiment with new visualization methods on the outputs of these analyses.
* To attempt to extract actionable information for public health responses to the ongoing pandemic.

## Current workflow
1. Sequences are bulk downloaded from the GISAID database.  All developers have signed the GISAID data access agreement, and sequences are not being re-distributed.
2. Sequences are aligned pairwise against the SARS-COV-2 reference genome using the Procrustean method implemented in [gotoh2](http://github.com/ArtPoon/gotoh2).  This module provides a method that progressively updates an existing alignment file with new sequence records, avoiding the re-alignment of previously released genomes.
3. Sequences are filtered for entries that are derived from non-human sources, and genomes that contain >5% fully ambiguous base calls (`N`s).
4. A pairwise genetic distance matrix is generated using [TN93](http://github.com/veg/tn93).
5. Clusters of effectively identical genome sequences are identified using the Python script `clustering.py`.  Representative sequences are written to a new FASTA file.  
6. A maximum likelihood phylogeny is reconstructed from the clustered FASTA file using the program [fasttree2](http://www.microbesonline.org/fasttree/) compiled for double precision.
