# limorhyde
`limorhyde` (linear models for rhythmicity, design) is an R package that enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [preprint](https://doi.org/10.1101/283622) and the [accompanying results](https://figshare.com/s/31dcb1346ef7f4268aa6).

## Install using drat
```R
install.packages('drat')
drat::addRepo('hugheylab')
install.packages('limorhyde', type = 'source')
```
You can update the package with `drat::addRepo('hugheylab')`, then `update.packages()`.

## Install using devtools
```R
install.packages('devtools')
devtools::install_github('hugheylab/limorhyde', build_vignettes = TRUE)
```
You can update the package using just the second line.

## Install using docker
You can use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed:
```
docker pull hugheylab/hugheyverse
```

## Getting started
See the documentation for the one function and the vignette:
```R
library('limorhyde')
?limorhyde
vignette(package = 'limorhyde')
```
