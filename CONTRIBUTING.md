
## Language translations

If you would like to contribute to CoVizu by providing a page translation in another language, thank you! 

* To get started, make a copy of `index.html` and add a [two-letter language code](https://en.wikipedia.org/wiki/ISO_639-1) suffix to the basename, *e.g.*, `index-fr.html`.

* We have moved all text elements into a single file `index.html` so that contributors do not have to search through the various JavaScript files for content to translate.
In order for this to work, we have had to create a [JSON](https://en.wikipedia.org/wiki/JSON) object at the top of the `<head>` section of `index.html`.
This maps keywords, *e.g.*, `number_cases`, to the language-specific replacement, *e.g.*, "nombre de cas".
A translation should be provided for every keyword in this JSON.

* Please check through the following HTML elements for content to translate:
  * `<label>` elements
  * `<option>` elements within `<select>` inputs
  * `title` tags within `<div>` and `<span>` elements
  * text enclosed in `<div>` and `<span>` elements
  * header elements, *e.g.*, `<h3>`

* There are several help texts to translate that are enclosed in `<p>` tags in the bottom half of the `index.html` document.
  * Text enclosed in the `title` tag of a `<span>` element must be written in a single line (no line breaks) to render properly.
  * Please leave the GISAID Acknowledgements statement: "We gratefully acknowledge all the Authors" in its original English.


## Development environment

Development of CoVizu is primarily being carried out on workstations and servers running Ubuntu 16.04+.
However, we have also run the system on desktop computers running macOS 10.13.
Most web development is tested in Google Chrome.

To setup the development environment:
1. Navigate to the directory containing `package.json`
2. Run the command `npm install`
3. Run the `run-server.sh` script


## Coding style

### Python

* Use four spaces to indent, and other conventions described in the [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide.  It is easier to use a Python-aware [IDE](https://en.wikipedia.org/wiki/Integrated_development_environment) like [PyCharm](https://www.jetbrains.com/pycharm/) or [Atom](https://atom.io/) to automate this for you. 

* Organize your code into functions to facilitate testing and so that methods can be called from other scripts.

* Isolate the main loop of your code under `if __name__ == '__main__'` to be executed only if the script is being run at the top level (from the command line). 

* Try to use `argparse` to provide help documentation for command line execution.

* Use lowercase with underscores to separate words for function and variable names.  Do not use [camel-case](https://en.wikipedia.org/wiki/Camel_case).

* Every function should open with a [docstring](https://www.python.org/dev/peps/pep-0257).  If the function is very brief with a small number of self-explanatory arguments, then a one-line docstring is fine.  Otherwise, use a multi-line docstring.  Use `:param varname: description` entries to document arguments. Use `:return: type, description` entries to document return values.   


## Committing your code

* Do not push commits to the `master` branch.  All development should be tracked in the `dev` branch and then merged into `master` after it has been tested.  
* If you are trying to implement something new that can potentially break parts of the code that other developers may be working on, then [create a new branch](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) and merge it into `dev` after you have finished building and testing the code. 
* Please try to write concise and informative [commit messages](https://xkcd.com/1296/).


## Front-end Testing 

* Front-end testing is carried out using [Cypress](https://www.cypress.io/). 
* For the purposes of this project, Cypress needs to be installed via direct download. An installation guide can be found [here](https://docs.cypress.io/guides/getting-started/installing-cypress#Direct-download). 
* The `baseUrl` is currently set as `"http://127.0.0.1:8001/`. Changes to it can be made in `Cypress.json`.

## Semantic Versioning

We follow [semantic versioning specifications](https://semver.org/). 

The version number (`MAJOR.MINOR.PATCH`) consist of the:
* `MAJOR` version - For changes that break backward compatibility. For example, modifying the JSON data structure required for running the front-end of CoVizu
* `MINOR` version - For backwards compatible features additions and other minor changes such as modifying documentation
* `PATCH` version - For backwards compatible bug fixes 

When opening a pull request, please specify the type of change (`major`, `minor` or `patch`) in the title. _e.g_ `patch: bug fixes for the search interface`.