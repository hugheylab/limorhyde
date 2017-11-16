# limorhyde

`limorhyde` (linear models for rhythmicity, design) is a package that enables differential rhythmicity analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we benchmarked it using simulations and real data from mice and humans, check out the [preprint]() and the [accompanying results]().

## Installation
Install the package in your local version of R:
```R
install.packages('devtools')
devtools::install_github('hugheylab/limorhyde')
```

Or use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies already installed:
```
docker pull hugheylab/hugheyverse
```

## Getting started
See the example in the documentation:
```R
library('limorhyde')
?limorhyde
```
