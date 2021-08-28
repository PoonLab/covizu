# openCoVizu: Open access real-time visualization of SARS-COV-2 genomic diversity

[CoVizu](https://github.com/PoonLab/covizu) is an ongoing open-source project to provide a public interface to visualize the global diversity of SARS-CoV-2 genomes in near real time.
Our specific objectives are (1) to process and visualize as much publicly available data as possible (*i.e.*, millions of genomes); (2) to reconstruct robust evolutionary and epidemiological relationships among these genomes; (3) to continually update outputs with new genomic data as frequently as possible, and; (4) to present this information in a rich and intuitive visual interface.

This `opendata` branch is a version of CoVizu that enables users to run the entire analysis and visualize results without the requirement of having special access to any particular database.
For example, you may want to try these methods on your own laboratory's data.
Hence, we provide script (`convert.py`) to adapt any combination of [FASTA](https://en.wikipedia.org/wiki/FASTA_format) and [tabular](https://en.wikipedia.org/wiki/Table_(information)) metadata into a common input format, and a second script to automate the analytical pipeline (`process.py`).

## Requirements
* [Python](https://www.python.org/) 3.6 or higher, and the following modules:
  * [BioPython](https://biopython.org/) version 1.7+
  * [mpi4py](https://pypi.org/project/mpi4py/)
  * [SciPy](https://www.scipy.org/) version 1.5+
* [minimap2](https://github.com/lh3/minimap2) version 2.1+ 
* [FastTree2](http://www.microbesonline.org/fasttree/) version 2.1.10+, compiled for [double precision](http://www.microbesonline.org/fasttree/#BranchLen)
* [TreeTime](https://github.com/neherlab/treetime) version 0.7.5+
* [RapidNJ](https://birc.au.dk/software/rapidnj/)
* [Pangolin](https://github.com/cov-lineages/pangolin/) (if your metadata file does not already include a field with Pango lineage classifications).

## Example usage

This sequence of commands demonstrates how to run the CoVizu pipeline using the [Open Data files](https://nextstrain.org/blog/2021-07-08-ncov-open-announcement) provided by the [Nextstrain](https://nextstrain.org/) development team, which are derived from the NCBI Genbank database.

1. First, we obtain the FASTA and metadata files from the Nextstrain data server:
   ```console
   $ wget https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz
   $ wget https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz
   ```
   Note, if you do not have the [`wget`](https://www.gnu.org/software/wget/) program then you can use a web browser to manually download the files at their respective URLs.

2. Next, we use the Python script `convert.py` to combine and convert these data files into a single xz-compressed JSON file.
   This script is designed to accommodate different tabular file formats for the metadata (*i.e.*, comma-, tab-separated) and different amounts of data.
   ```console
   $ python3 convert.py sequences.fasta.xz metadata.tsv.gz --xz --mgz --region region --division division --outfile opendata.json.xz
   ```
   The `--xz` and `--mgz` flags tell the script that the FASTA and metadata files are xz- and gzip-compressed, respectively.
   Leaving these files in their compressed state while streaming data minimizes data storage requirements.
   Converting these files (about 1.2 million genome records) consumed about an hour on my AMD/Ryzen workstation running Ubuntu.

3. To run the analysis, we call the Python script `process.py` on the JSON file produced by `convert.py`:
   ```console
   $ python3 process.py opendata.2021-08-25.xz --ft2bin fasttree2-mp -np 8
   ```
   The `-np` argument is used to set the number of cores for [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) processing.
   You should change this number according to your hardware (number of cores, RAM).
   In this example, I have used the `ft2bin` flag to specify the non-standard name of the [Fasttree2](http://www.microbesonline.org/fasttree/) binary on my computer.
   
   For larger databases, you may need to launch the pipeline in the background and reduce the number of cores:
   ```console
   nohup python3 process.py opendata.2021-08-25.xz -np 4 & disown
   ```
   Running on four cores consumed about 18 Gb of RAM and took about 14 hours to complete for this database of 1.2 million genomes.

4. The analysis script in step 3 writes three time-stamped output files to sub-directory `data/`.  To view these data in the CoVizu web interface, you need to copy these files as follows:
   ```console
   $ cp data/clusters.2021-08-25T17:24:55.json data/clusters.json
   $ cp data/dbstats.2021-08-25T17:24:55.json data/dbstats.json
   $ cp data/timetree.2021-08-25T17:24:55.nwk data/timetree.nwk
   ```
   Keeping the time-stamped data files is useful for revisiting results at a previous time, and they are relatively compact.

5. Launch a local webserver with `bash run-server.sh` or `python3 -m http.server 8001 --bind 127.0.0.1`, and navigate your browser to `localhost:8001`.


## Acknowledgements
The development and validation of these scripts was made possible by the labs who have generated and contributed SARS-COV-2 genomic sequence data, much of which has been curated and published by the [GISAID Initiative](https://www.gisaid.org/).  We sincerely thank these labs for making this information available to the public and open science.
