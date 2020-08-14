
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
0 0 * * * nohup python /home/covid/autobot.py >> /home/covid/Autobot.log 2>&1
```
Note that the absolute paths in the above are specific to the filesystem on our own server, and that you should replace this with your own path specification!


## Coding style

### Python

* Use four spaces to indent, and other conventions described in the [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide.  It is easier to use a Python-aware [IDE](https://en.wikipedia.org/wiki/Integrated_development_environment) like [PyCharm](https://www.jetbrains.com/pycharm/) or [Atom](https://atom.io/) to automate this for you. 

* Organize your code into functions to facilitate testing and so that methods can be called from other scripts.

* Isolate the main loop of your code under `if __name__ == '__main__'` to be executed only if the script is being run at the top level (from the command line). 

* Try to use `argparse` to provide help documentation for command line execution.

* Use lowercase with underscores to separate words for function and variable names.  Do not use [camel-case](https://en.wikipedia.org/wiki/Camel_case).

* Every function should open with a [docstring](https://www.python.org/dev/peps/pep-0257).  If the function is very brief with a small number of self-explanatory arguments, then a one-line docstring is fine.  Otherwise, use a multi-line docstring.  Use `:param varname: description` entries to document arguments. Use `:return: type, description` entries to document return values.   


### R 

* Indent with two spaces
* Use base R whenever possible
* Use `.` to separate words in variable and functino names, not `_`
* Use `#'` prefix to document functions, *e.g.*:
  ```R
  #' @param node: str, label of current node variant
  #' @param parent: str, label of current node's parental variant
  #' @param el: str, edge list from minimum spanning tree
  #' @return linearized vector of parent->child pairs
  traverse <- function(node, parent, el, edges=c()) {
  ```
* Place any package requirements (*i.e.*, `require(igraph)`) at the top of the script


## Commiting your code

Do not push commits to the `master` branch.  All development should be tracked in the `dev` branch and then merged into `master` after it has been tested.  

If you are trying to implement something new that can potentially break parts of the code that other developers may be working on, then [create a new branch](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) and merge it into `dev` after you have finished building and testing the code. 

Please try to write concise and informative [commit messages](https://xkcd.com/1296/).

