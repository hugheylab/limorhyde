# limorhyde
`limorhyde` (linear models for rhythmicity, design) is an R package that enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [preprint]() and the [accompanying results]().

## Install using drat
```R
install.packages('drat')
drat::addRepo('hugheylab')
install.packages('limorhyde', type='source')
```
You can update the package with `drat::addRepo('hugheylab')`, then `update.packages()`.

## Install using devtools
```R
install.packages('devtools')
devtools::install_github('hugheylab/limorhyde')
```
You can update the package using just the second line.

## Install using docker
You can use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed:
```
docker pull hugheylab/hugheyverse
```

## Getting started
See the example in the documentation:
```R
library('limorhyde')
?limorhyde
```
