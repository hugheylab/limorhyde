## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## ---- message = FALSE----------------------------------------------------
library(annotate)
library(limorhyde)
library(GEOquery)
library(dplyr)
library(tidyr)
library(knitr)
library(limma)
library(rain)
library(ggplot2)
library(org.Mm.eg.db)

source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))

## ------------------------------------------------------------------------
period = 24
qvalRhyCutoff = 0.15
qvalDrCutoff = 0.1

