## Data

If you want to contribute to this project, please create an account with GISAID so you sign their data access agreement.
GISAID sequence data was retrieved by Emmanuel on April 3 and is being stored on Langley at `/home/covid/gisaid_cov2020_sequences.fasta`.
Currently, there are 129 genome sequences that were collected in Canada.

## Analysis plan

### Alignment

Generally we need a multiple sequence alignment if we want to use most conventional methods for comparative analysis or phylogenetic reconstruction.
Since there are roughly 4,000 sequences of about 30,000 nt each, building an MSA will be a very time-consuming process.
I think our first problem is to develop a fast pairwise-based method for generating an MSA.  
See the `cutter.py` script in the PoonLab repo.

### Phylogeny

With an MSA, it will be straight-forward to generate a tree by maximum likelihood.  A lot of identical sequences will be collapsed.  If we use IQ-TREE (currently my preferred method), we should disable the model selection phase and just use a simple TN93 model (or check the literature to see what named model is supported).

### Annotation

In order to reconstruct transmission events between countries, we need to annotate the phylogeny with these metadata.  Sample collection information is embedded in the sequence labels.  

