# CoVizu: Real-time visualization of SARS-COV-2 genomic diversity

CoVizu is an open source project to develop a `near real time' SARS-CoV-2 genome analysis and visualization system that highlights potential cases of importation from other countries or ongoing community transmission.

The current mode of visualization employed by CoVizu that we are tentatively referring to as a "beadplot":
![](vignettes/beadplot.png)


### How to read a beadplot:

* Each horizontal line segment represents a unique SARS-CoV-2 genomic sequence variant.  The emergence of a single new mutation in an infection is sufficient to establish a new variant.  A given variant may be observed multiple times as identical genome sequences, where `identity' is loosely defined to accommodate genomes with incomplete coverage and ambiguous base calls.  (Following GISAID's definition of a "complete genome", we require a minimum sequence length of 29,000 nt.)

* Each circle represents one or more cases of a variant that were sampled on a given date.  The size (area) of the circle is proportional to the number of sequences.

* Cases are arranged in chronological order from left to right.

* Circles are coloured red if *any* of the genomes was sampled in Canada on that date.

  > The first red circle in a series of cases on the same horizontal line represents the first sample of that variant in Canada and implies a new importation.  A series of red circles on the same line implies community transmission.

* Variants are labelled with the name of the earliest genome record.  To reduce clutter, however, it may become necessary to filter labels to only variants that have a minimum number of cases.

* Vertical lines connect variants that are related by a minimum spanning tree, which gives a *rough* approximation of transmission links.  The variant at the bottom terminus of the vertical line is the putative source.  

**It is not feasible to reconstruct accurate links using only genomic data.**  However, our objective is to identify population-level events like importations into Canada, not to attribute a transmission to a specific source individual.

* Circles that are the first "bead" on the horizontal line are related to the same ancestral variant if they intersect the same vertical line.  This scenario implies that multiple lineages descend from the same ancestor.

  > If none of the ancestral cases were sampled in Canada, then a first "bead" of a new variant that is red implies importation and mutation.

* The relative location of variants along the vertical axis does not convey any information.  The variants are sorted with respect to the vertical axis such that ancestral variants are always below their "descendant" variants.


## Current workflow

1. Sequences are bulk downloaded from the GISAID database.  All developers have signed the GISAID data access agreement, and sequences are not being re-distributed.

2. Sequences are aligned pairwise against the SARS-COV-2 reference genome using the Procrustean method implemented in [gotoh2](http://github.com/ArtPoon/gotoh2) - see `updater.py`.  This module provides a method that progressively updates an existing alignment file with new sequence records, avoiding the re-alignment of previously released genomes.

3. Sequences are filtered using `filtering.py` for entries that are derived from non-human sources, incomplete genomes, and genomes that contain >5% fully ambiguous base calls (`N`s).

4. A pairwise genetic distance matrix is generated using [TN93](http://github.com/veg/tn93) - only distances below a cutoff of `0.0001` are recorded to the output file.

5. Identical sequences are collapsed into a small number of variants using the Python script `variants.py`.  Note that this does not require that sequences overlap end-to-end, nor that the overlapping regions are an exact match (accommodating ambiguous base calls, for example).

6. A pairwise TN93 genetic distance matrix is generated for the reduced data set of unique genomic variants.

7. Variants are clustered with respect to the pairwise distance matrix using hierarchical clustering.  A [minimum spanning tree](https://en.wikipedia.org/wiki/Minimum_spanning_tree) is reconstructed for each cluster, and rooted at the variant that is closest to the earliest sampled variant (Wuhan, IPBCAMS-WH-01).  The results are written to [JSON](https://en.wikipedia.org/wiki/JSON) files.

8. Visualizations are genearted from the JSON data using the [d3.js](https://en.wikipedia.org/wiki/D3.js) Javascript framework.


## Acknowledgements
The development and validation of these scripts was made possible by the labs who have generated and contributed SARS-COV-2 genomic sequence data that is curated and published by [GISAID](https://www.gisaid.org/).  We sincerely thank these labs for making this information available to the public and open science.
