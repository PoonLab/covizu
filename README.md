# CoVizu: Open access real-time visualization of SARS-COV-2 genomic diversity
## VirusSeq version

[CoVizu](https://github.com/PoonLab/covizu) is an ongoing open-source project to provide a public interface to visualize the global diversity of SARS-CoV-2 genomes in near real time.
Our specific objectives are (1) to process and visualize as much publicly available data as possible (*i.e.*, millions of genomes); (2) to reconstruct robust evolutionary and epidemiological relationships among these genomes; (3) to continually update outputs with new genomic data as frequently as possible, and; (4) to present this information in a rich and intuitive visual interface.


## Requirements
* [Python](https://www.python.org/) 3.6 or higher, and the following modules:
  * [BioPython](https://biopython.org/) version 1.7+
  * [mpi4py](https://pypi.org/project/mpi4py/)
  * [SciPy](https://www.scipy.org/) version 1.5+
* [minimap2](https://github.com/lh3/minimap2) version 2.1+ 
* [FastTree2](http://www.microbesonline.org/fasttree/) version 2.1.10+, compiled for [double precision](http://www.microbesonline.org/fasttree/#BranchLen)
* [TreeTime](https://github.com/neherlab/treetime) version 0.7.5+
* [RapidNJ](https://birc.au.dk/software/rapidnj/)

## Example usage

This sequence of commands demonstrates how to run the CoVizu pipeline using:
* [Open Data files](https://nextstrain.org/blog/2021-07-08-ncov-open-announcement) provided by the [Nextstrain](https://nextstrain.org/) development team, which are derived from the NCBI Genbank database.
* [Canadian VirusSeq Data Files](https://virusseq-dataportal.ca/)

1. First, we obtain the FASTA and metadata files from the Nextstrain data server:
   ```console
   $ wget https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz
   $ wget https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
   ```
   Note, if you do not have the [`wget`](https://www.gnu.org/software/wget/) program then you can use a web browser to manually download the files at their respective URLs.

2. Next, we download the FASTA and metadata files from the VirusSeq portal database and the PANGO lineage classifications from the [Viral AI](https://viral.ai/collections) database:
   ```console
   $ python3 retrieve.py data
   ``` 
3. To run the analysis, we call the Python script `process.py`:
   ```console
   $ python3 process.py -np 8 --vsfasta data/virusseq.fasta.xz --vsmeta data/virusseq.metadata.tsv.gz --vspango data/viralai.csv --opfasta data/sequences.fasta.xz --opmeta data/metadata.tsv.gz --ft2bin fasttree2-mp
   ```
   The `-np` argument is used to set the number of cores for [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) processing.
   You should change this number according to your hardware (number of cores, RAM).
   `--vsfasta` and `--vsmeta` are used to provide the paths to the FASTA and metadata files from the VirusSeq portal. 
   `--vspango` is used to provide the path to PANGO lineage classifications from the Viral AI database
   `--opfasta` and `--opmeta` are used to provide the paths to the FASTA and metadata files from the Nextstrain data server.
   In this example, I have used the `ft2bin` flag to specify the non-standard name of the [Fasttree2](http://www.microbesonline.org/fasttree/) binary on my computer.

4. The analysis script in step 3 writes three time-stamped output files to sub-directory `data/`.  To view these data in the CoVizu web interface, you need to copy these files as follows:
   ```console
   $ cp data/clusters.2021-08-25T17:24:55.json data/clusters.json
   $ cp data/dbstats.2021-08-25T17:24:55.json data/dbstats.json
   $ cp data/timetree.2021-08-25T17:24:55.nwk data/timetree.nwk
   ```
   Keeping the time-stamped data files is useful for revisiting results at a previous time, and they are relatively compact.

5. Launch the Node.js Web Server

   If you're running the server for the first time, navigate to the directory containing `package.json` and run the command `npm install` to install all the required dependencies to run the server.

   Launch the local webserver with `npm start`, allow up to a minute for the server to initialize and navigate your browser to `localhost:8001`.


## Acknowledgements
The development and validation of these scripts was made possible by the labs who have generated and contributed SARS-COV-2 genomic sequence data, much of which has been curated and published by the [GISAID Initiative](https://www.gisaid.org/).  We sincerely thank these labs for making this information available to the public and open science.

Computational support provided by the [McArthur Lab](https://mcarthurbioinformatics.ca/) (Michael G. DeGroote Institute for Infectious Disease Research, McMaster University).

PANGO lineage classifications were retrieved from the [Viral AI](https://viral.ai/collections) database.
