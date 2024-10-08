

# mglm4twin 0.5.0

The `mglm4twin` package fits multivariate generalized linear models 
for twin and family data. (Bonat and Hjelmborg, 2020). Check the paper
[here](https://link.springer.com/article/10.1007/s10519-021-10095-3).


## Introduction

`mglm4twin` fits multivariate generalized linear models for twin and 
family data. It allows use a different linear predictor for each 
response variable of a multivariate response. 
The response variable can be continuous or discrete, like counts and 
binomial and also limited continuous or discrete/continuous inflated 
responses. The most important and relevant feature is that many covariance 
structures can be used to model the relations among different traits.

## Download and install

### Linux/Mac

Use the `devtools` package (available from
[CRAN](http://cran-r.c3sl.ufpr.br/web/packages/devtools/index.html)) to
install automatically from this GitHub repository:


``` r
library(devtools)
install_github("bonatwagner/mglm4twin")
```

Alternatively, download the package tarball: [mglm4twin_0.5.0.tar.gz][]
and run from a UNIX terminal (make sure you are on the container file
directory):


```
R CMD INSTALL -l /path/to/your/R/library mglm4twin_0.5.0.tar.gz
```

Or, inside an `R` session:


```
install.packages("mglm4twin_0.5.0.tar.gz", repos = NULL,
                 lib.loc = "/path/to/your/R/library",
                 dependencies = TRUE)
```

Note that `-l /path/to/your/R/library` in the former and `lib.loc =
"/path/to/your/R/library"` in the latter are optional. Only use it if
you want to install in a personal library, other than the standard R
library.

### Windows

Download Windows binary version: [mglm4twin_0.5.0.zip][] (**do not unzip
it under Windows**), put the file in your working directory, and from
inside `R`:


```
install.packages("mglm4twin_0.5.0.zip", repos = NULL,
                 dependencies = TRUE)
```


## Documentation

The reference manual in PDF can be found here: [mglm4twin-manual.pdf][]

## Contributing

This R package is develop using [`roxygen2`][] for documentation and
[`devtools`] to check and build. Also, we adopt the [Gitflow worflow][]
in this repository. Please, see the
[instructions for contributing](./CONTRIBUTING.md) to collaborate.

## License

This package is released under the
[GNU General Public License (GPL) v3.0][].

See [LICENSE](./LICENSE)

<!-- links -->



[GNU General Public License (GPL) v3.0]: http://www.gnu.org/licenses/gpl-3.0.html
[`roxygen2`]: https://github.com/klutometis/roxygen
[`devtools`]: https://github.com/hadley/devtools
[mglm4twin_0.5.0.tar.gz]: http://www.leg.ufpr.br/~wagner/mglm4twin/source/mglm4twin_0.5.0.tar.gz
[mglm4twin_0.5.0.zip]: http://www.leg.ufpr.br/~wagner/mglm4twin/source/mglm4twin_0.5.0.zip
[mglm4twin-manual.pdf]: http://www.leg.ufpr.br/~wagner/mglm4twin/source/mglm4twin-manual.pdf
[Gitflow workflow]: http://nvie.com/posts/a-successful-git-branching-model/
[Wagner Hugo Bonat]: http://www.leg.ufpr.br/~wagner
