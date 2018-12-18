# limorhyde
`limorhyde` (linear models for rhythmicity, design) is an R package that enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [paper](https://doi.org/10.1177/0748730418813785) and the [accompanying results](https://doi.org/10.6084/m9.figshare.5945569).

## Installation
First install drat.
```R
install.packages('drat')
```

Then add the following line to your `.Rprofile` file (located at "~/.Rprofile"), which gets run every time R starts. See [here](https://csgillespie.github.io/efficientR/3-3-r-startup.html#r-startup) for details.
```R
drat::addRepo('hugheylab')
```

Now you can install the package.
```R
install.packages('limorhyde', type = 'source')
```
You can update the package using `update.packages()`.

## Docker
You can also use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.
```bash
docker pull hugheylab/hugheyverse
```

## Getting started
See the vignette and the documentation.
```R
browseVignettes('limorhyde')
?limorhyde::limorhyde
```
