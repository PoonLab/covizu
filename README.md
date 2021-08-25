# CoVizu: Real-time visualization of SARS-COV-2 genomic diversity

CoVizu is an open source project to develop a `near real time' SARS-CoV-2 genome analysis and visualization system that highlights potential cases of importation from other countries or ongoing community transmission.

This `opendata` branch is a version of CoVizu that enables users to run the entire analysis and visualize results without the requirement of special access to any particular database.  

## Requirements
* [Python](https://www.python.org/) 3.6 or higher, and the following modules:
  * [BioPython](https://biopython.org/) version 1.7+
  * [mpi4py](https://pypi.org/project/mpi4py/)
  * [SciPy](https://www.scipy.org/) version 1.5+
* [minimap2](https://github.com/lh3/minimap2) version 2.1+ 
* [FastTree2](http://www.microbesonline.org/fasttree/) version 2.1.10+, compiled for [double precision](http://www.microbesonline.org/fasttree/#BranchLen)
* [TreeTime](https://github.com/neherlab/treetime) version 0.7.5+
* [RapidNJ](https://birc.au.dk/software/rapidnj/)
* [Pangolin](https://github.com/cov-lineages/pangolin/)

## Usage

This sequence of commands demonstrates how to run the CoVizu pipeline using the [Open Data files](https://nextstrain.org/blog/2021-07-08-ncov-open-announcement) provided by the Nextstrain development team, which are derived from the NCBI Genbank database.

1. First, we obtain the FASTA and metadata files from the Nextstrain data server:
   ```console
   $ wget https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz
   $ wget https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
   ```
   Note, if you do not have the [`wget`](https://www.gnu.org/software/wget/) program then you can use a web browser to manually download the files at their respective URLs.

2. Next, we use the Python script `convert.py` to combine and convert these data files into a single xz-compressed JSON file.
   This script is designed to accommodate different tabular file formats for the metadata, *i.e.*, comma-separated, tab-separated, and different amounts of data.
   ```console
   $ python3 convert.py sequences.fasta.xz metadata.tsv.gz --xz --mgz --region region --division division --outfile opendata.json.xz
   ```
   The `--xz` and `--mgz` flags tell the script that the FASTA and metadata files are xz- and gzip-compressed, respectively.
   Leaving these files in their compressed state while streaming data minimizes data storage requirements.
   Converting these files (about 1.2 million genome records) consumed about 45 minutes on my Ubuntu workstation.

3. To run the analysis, we call the Python script `process.py` on the JSON file produced by `convert.py`:
   ```console
   $ python3 process.py opendata.json.xz
   ```


## Acknowledgements
The development and validation of these scripts was made possible by the labs who have generated and contributed SARS-COV-2 genomic sequence data that is curated and published by [GISAID](https://www.gisaid.org/).  We sincerely thank these labs for making this information available to the public and open science.
