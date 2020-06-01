
# Development environment

Development of CoVizu is primarily being carried out on workstations and servers running Ubuntu 16.04+. 
However, we have also run the system on desktop computers running macOS 10.13.

## Dependencies

* a C build environment
* Mozilla [geckodriver](https://github.com/mozilla/geckodriver) v0.26+
* [Python](https://www.python.org/) 3.6 or higher, and the following modules:
  * [Selenium](https://github.com/SeleniumHQ/selenium/) version 3.14.1+ 
  * [gotoh2](https://github.com/ArtPoon/gotoh2/)
  * [networkx](https://networkx.github.io/) version 2.3+
  * [BioPython](https://biopython.org/) version 1.7+
* GNU [sed](https://www.gnu.org/software/sed/) stream editor
* [TN93](https://github.com/veg/tn93) v1.0.6
* [FastTree2](http://www.microbesonline.org/fasttree/) version 2.1.10+, compiled for [double precision](http://www.microbesonline.org/fasttree/#BranchLen)
* [TreeTime](https://github.com/neherlab/treetime) version 0.7.5+
* [R](https://cran.r-project.org/) 3.6+, and the following packages:
  * [igraph](https://igraph.org/r/) version 1.2+
  * [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html) version 1.6+
  * [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html) version 0.15


## GISAID database access

To retrieve data from GISAID, you must agree to the terms of their [Database Access Agreement](https://www.gisaid.org/registration/terms-of-use/) and then complete and submit the [registration form](https://www.gisaid.org/registration/register/).
Note that access credentials are only necessary if you are doing any back-end development. 
For front-end development, we upload the necessary JSON files to our GitHub repository.
These files do not contain any genome sequence data - only sequence labels and associated metadata (date and country of sample collection), and results from our clustering and phylogenetic reconstruction analyses are stored.

The following environment variables need to be defined for the automated GISAID download scripts (`SeliniumAutobot.py` and ), which can be done by adding the following lines to the `.bashrc` shell script in your home directory, for example: 
```bash
export gisaid_u_variable='<your GISAID username>'
export gisaid_pw_variable='<your GISAID password>'
```

This does mean that these database access credentials are written to a plain text file in your filesystem.
Consequently, you should only use this approach if you are running CoVizu in a fairly secure computing environment (*e.g.*, a physically- and password-secured workstation, ideally with multi-factor authentication). 


The downloading scripts can be automated through `crontab` on Linux:
```
0 0 * * * nohup python /home/covid/SeliniumAutobot.py >> /home/covid/Autobot.log 2>&1
```


# Workflow

## Back-end

### Data collection

Presently, we are using the `SeleniumAutobot.py` script to retrieve any records in the GISAID database with an upload date in the last 24 hours.
We limit the number of records to reduce the burdern on the GISAID database servers and to minimize sample processing time.
Genome records that were retrieved in a previous query do not have to be re-aligned.

### Sequence alignent



### Data processing

## Front-end

### Data serialization



