# limorhyde

[![check-coverage-deploy](https://github.com/hugheylab/limorhyde/workflows/check-coverage-deploy/badge.svg)](https://github.com/hugheylab/limorhyde/actions)
[![codecov](https://codecov.io/gh/hugheylab/limorhyde/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/limorhyde)

`limorhyde` (linear models for rhythmicity, design) enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [paper](https://doi.org/10.1177/0748730418813785) and the [accompanying results](https://doi.org/10.6084/m9.figshare.5945569).

## Installation

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter:

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once. Then you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('limorhyde')
```

Alternatively, you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('limorhyde', site_repository = 'https://hugheylab.github.io/drat/')
```

## Usage

For an introduction to the package, read the [vignette](https://limorhyde.hugheylab.org/articles/introduction.html). For more details, check out the [reference documentation](https://limorhyde.hugheylab.org/reference/index.html).
