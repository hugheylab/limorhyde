# limorhyde
`limorhyde` (linear models for rhythmicity, design) is an R package that enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [preprint](https://doi.org/10.1101/283622) and the [accompanying results](https://figshare.com/s/31dcb1346ef7f4268aa6).

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
See the documentation and the vignette.
```R
library('limorhyde')
?limorhyde
vignette('introduction', package = 'limorhyde')
```
