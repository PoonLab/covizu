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
* [Node.js](https://nodejs.org/en/download/)
* [npm](https://docs.npmjs.com/about-npm-versions)

If running locally (without dedicated GISAID feed):
* [Pangolin](https://github.com/cov-lineages/pangolin/)

## Initialize Development Environment
If you're running the server for the first time run the following commands:
```
 git clone https://github.com/PoonLab/covizu
 cp timetree.nwk covizu/data/
 cp dbstats.json covizu/data/
 cp clusters.json covizu/data/
 cd covizu
 npm install
 npm start
```
Once you launch the local webserver with `npm start`, allow up to a minute for the server to initialize and then navigate your browser to localhost:8001