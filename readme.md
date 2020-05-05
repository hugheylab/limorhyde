# limorhyde

`limorhyde` (linear models for rhythmicity, design) enables differential analysis of circadian transcriptome data using state-of-the-art methods for differential expression.

For details about the method and to see how we used it to analyze circadian transcriptome data from various experimental designs, check out the [paper](https://doi.org/10.1177/0748730418813785) and the [accompanying results](https://doi.org/10.6084/m9.figshare.5945569).

## Installation

First add the hugheylab repository to your repos. There are multiple ways to do this.

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter the following values.

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once.

Alternatively, enter the following command each time you want to install or update the package.

```R
options(repos = c(getOption('repos'), 'https://hugheylab.github.io/drat/'))
```

Now you can install the package.

```R
install.packages('limorhyde', type = 'source')
```

You can update the package using `update.packages()`.

There's also [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.

```bash
docker pull hugheylab/hugheyverse
```

## Usage

For an introduction to the package, read the [vignette](articles/introduction.html). For more details, check out the [reference documentation](reference/index.html).
