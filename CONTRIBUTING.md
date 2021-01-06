
# Development environment

Development of CoVizu is primarily being carried out on workstations and servers running Ubuntu 16.04+. 
However, we have also run the system on desktop computers running macOS 10.13.


## Coding style

### Python

* Use four spaces to indent, and other conventions described in the [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide.  It is easier to use a Python-aware [IDE](https://en.wikipedia.org/wiki/Integrated_development_environment) like [PyCharm](https://www.jetbrains.com/pycharm/) or [Atom](https://atom.io/) to automate this for you. 

* Organize your code into functions to facilitate testing and so that methods can be called from other scripts.

* Isolate the main loop of your code under `if __name__ == '__main__'` to be executed only if the script is being run at the top level (from the command line). 

* Try to use `argparse` to provide help documentation for command line execution.

* Use lowercase with underscores to separate words for function and variable names.  Do not use [camel-case](https://en.wikipedia.org/wiki/Camel_case).

* Every function should open with a [docstring](https://www.python.org/dev/peps/pep-0257).  If the function is very brief with a small number of self-explanatory arguments, then a one-line docstring is fine.  Otherwise, use a multi-line docstring.  Use `:param varname: description` entries to document arguments. Use `:return: type, description` entries to document return values.   


## Committing your code

Do not push commits to the `master` branch.  All development should be tracked in the `dev` branch and then merged into `master` after it has been tested.  

If you are trying to implement something new that can potentially break parts of the code that other developers may be working on, then [create a new branch](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) and merge it into `dev` after you have finished building and testing the code. 

Please try to write concise and informative [commit messages](https://xkcd.com/1296/).

