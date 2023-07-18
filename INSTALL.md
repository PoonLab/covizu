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
* [Node.js](https://nodejs.org/en/download/) version 18.0.0+
* [npm](https://docs.npmjs.com/about-npm-versions)

If running locally (without dedicated GISAID feed):
* [Pangolin](https://github.com/cov-lineages/pangolin/)


## Installation

* Navigate to the directory in your filesystem under which you want to install the `covizu` directory

* Clone the repository:
```
git clone https://github.com/PoonLab/covizu
```

* Install with:
```
sudo python3 setup.py install  # omit `sudo` if you don't have superuser privileges
```

If you're running the server (front-end) for the first time, obtain the following data files from our main server at 
https://filogeneti.ca/covizu/data/:
* `timetree.nwk`
* `dbstats.json`
* `clusters.json`
and save your local copies under `covizu/data/`.

Next, make 3 separate environment files for the NodeJS server. 
#### `.env.test`
```
HTTP_PORT='8001'
NODE_ENV='TEST'
DATA_FOLDER='data_test'
```
#### `.env.dev`
```
HTTP_PORT='8001'
NODE_ENV='DEV'
DATA_FOLDER='data'
ADMIN_USERNAME='admin'
ADMIN_PASSWORD='password'
COVIZU_USERNAME='covizu'
COVIZU_PASSWORD='password'
DB_URL='localhost:27017'
```
#### `.env.prod` (if running a production server with SSL enabled)
```
HTTP_PORT='8001'
HTTPS_PORT='8002'
NODE_ENV='PROD'
DATA_FOLDER='data'
PRVTKEY='/path/to/private-key-file.pem'
CRT='/path/to/certificate-file.crt'
ADMIN_USERNAME='admin'
ADMIN_PASSWORD='password'
COVIZU_USERNAME='covizu'
COVIZU_PASSWORD='password'
DB_URL='localhost:27017'
```

Install and configure the `MongoDB` database by following these [instructions](DBINSTALL.md)

Finally, run the following commands:
```
 cd covizu
 npm install
 npm run start-db1
```

Once you launch the local webserver with `npm start`, allow up to a minute for the server to initialize and then navigate your browser to `localhost:8001`.

To run the end-to-end tests with [Cypress](http://cypress.io) start the test server
```
npm test
```
and in a new terminal terminal window run 
```
npx cypress run
```

