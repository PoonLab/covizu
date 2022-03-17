## Dependencies

* [Python](https://www.python.org/) 3.6 or higher, and the following modules:
  * [BioPython](https://biopython.org/) version 1.7+
  * [mpi4py](https://pypi.org/project/mpi4py/)
  * [SciPy](https://www.scipy.org/) version 1.5+
* [minimap2](https://github.com/lh3/minimap2) version 2.1+ 
* [FastTree2](http://www.microbesonline.org/fasttree/) version 2.1.10+, compiled for [double precision](http://www.microbesonline.org/fasttree/#BranchLen)
* [TreeTime](https://github.com/neherlab/treetime) version 0.7.5+
* [RapidNJ](https://birc.au.dk/software/rapidnj/)
* [git](https://git-scm.com/)

If running locally (without dedicated GISAID feed):
* [Pangolin](https://github.com/cov-lineages/pangolin/)

## Installation

* Navigate to the directory in your filesystem under which you want to install the `covizu` directory
* Clone the repository:
```
git clone https://github.com/PoonLab/covizu
```
* Enter the new directory and install the submodules `pango-designation` and `ProblematicSites_SARS-CoV2`:
```
cd covizu
git submodule init
git submodule update
```
* Install the package into your Python module directory
```
sudo python3 setup.py install  # omit `sudo` if you don't have superuser privileges
```
