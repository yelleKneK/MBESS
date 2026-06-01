# MBESS

**MBESS** (pronounced 'em-bes') is an R package that implements methods
useful in designing research studies and analyzing data, with particular
emphasis on methods developed for or used within the behavioral,
educational, and social sciences (broadly defined). Many of the methods,
however, apply across a wide variety of disciplines.

MBESS was originally an acronym for 'Methods for the Behavioral,
Educational, and Social Sciences.' As the package became more general it
became an orphan acronym: what was once an acronym is now literally its
name.

## Capabilities

MBESS provides suites of functions for:

* effect sizes (standardized and unstandardized);
* confidence intervals for effect sizes, including noncentral methods;
* sample size planning from the accuracy in parameter estimation (AIPE),
  power analytic, equivalence, and minimum risk point estimation
  perspectives;
* mediation analysis;
* properties of distributions; and
* a variety of utility functions.

See the help pages (for example, `library(help = 'MBESS')`) for the full
list of functions.

## Installation

Install the released version from CRAN:

```r
install.packages('MBESS')
```

The development version is available on GitHub:

```r
# install.packages('remotes')
remotes::install_github('yelleKneK/MBESS')
```

## Citation

To cite MBESS in publications, run:

```r
citation('MBESS')
```

The primary references are Kelley (2007, *Behavior Research Methods*,
39, 979-984) and Kelley (2007, *Journal of Statistical Software*, 20(8),
1-24).

## Acknowledgments

MBESS has benefited from the contributions of many people. See the
`CONTRIBUTORS` file that ships with the package, located via
`system.file('CONTRIBUTORS.md', package = 'MBESS')`.

## Maintenance status

MBESS is mature, stable, and widely used, and it remains fully
supported. Maintenance mode does not mean the package is finished or
abandoned: MBESS will continue to be maintained so that it keeps
working with current versions of R and its dependencies, and it will
remain on CRAN. Existing code that relies on MBESS will keep working,
so users can continue to depend on it with confidence.

New methodological development now happens in its successor package,
**DMAR** (Design, Measurement, and Analysis in R), which is available
on CRAN. DMAR is where new methods will appear, and MBESS and DMAR are
designed to coexist, so there is no need to migrate existing work.

## License

GPL-2 | GPL-3. Maintainer: Ken Kelley (kkelley@nd.edu).
