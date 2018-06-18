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

## ---- message = FALSE----------------------------------------------------
gse = getGEO('GSE34018')

## ---- message = FALSE----------------------------------------------------
sm = as_tibble(pData(phenoData(gse[[1]]))) %>%
  transmute(title = title,
            sample = geo_accession,
            genotype = `genotype/variation:ch1`,
            time = as.numeric(sub("_.*", "", sub(".*_ZT", "", title))),
            cond = factor(genotype, c('wild-type', 'Reverb alpha/beta DKO'),
                          c('wild-type', 'knockout'))) %>%
  arrange(cond, time)

kable(sm[1:5,])

## ------------------------------------------------------------------------
sm = bind_cols(sm, limorhyde(sm, 'time', period))

## ----message = FALSE-----------------------------------------------------
mapping = getGeneMapping(featureData(gse$GSE34018_series_matrix.txt.gz))
emat = log2(calcExprsByGene(gse$GSE34018_series_matrix.txt.gz, mapping))[,sm$sample]

## ---- message = FALSE----------------------------------------------------
rhyRain = rainWrapper(sm, emat, period)

rhyGenes = rhyRain %>%
  group_by(geneId) %>%
  summarize(pVal = min(pVal)) %>%
  group_by() %>%
  mutate(qval = p.adjust(pVal, method = 'BH')) %>%
  arrange(qval)

kable(rhyGenes[1:5,])

## ---- message = FALSE----------------------------------------------------
design = model.matrix(~ cond * (time_cos + time_sin), data = sm)

fit = lmFit(emat, design)
fit = eBayes(fit, trend=TRUE)
drLimma = topTable(fit, coef = 5:6, number = Inf)

drLimma$geneId = rownames(drLimma)
rownames(drLimma) = NULL
drLimma = semi_join(drLimma, filter(rhyGenes, qval <= qvalRhyCutoff), by = 'geneId')
drLimma$adj.P.Val = p.adjust(drLimma$P.Value, method = 'BH')

kable(drLimma[1:5,])

## ---- message = FALSE----------------------------------------------------
design = model.matrix(~ cond + time_cos + time_sin, data = sm)

fit = lmFit(emat, design)
fit = eBayes(fit, trend = TRUE)
deLimma = topTable(fit, coef = 2, number = Inf)

deLimma$geneId = rownames(deLimma)
rownames(deLimma) = NULL
deLimma = anti_join(deLimma, filter(drLimma, adj.P.Val <= qvalDrCutoff), by = 'geneId')
deLimma$adj.P.Val = p.adjust(deLimma$P.Value, method = 'BH')

kable(deLimma[1:5,])

## ----fig.width = 4, fig.height = 3, message = FALSE----------------------
ggplot(data = deLimma) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), size = 0.2, alpha = 0.5) +
  labs(x = expression(log[2]*' fold-change'), y = expression(-log[10]*' '*q[DE]))

## ----fig.width = 5, fig.height = 4, message = FALSE----------------------
geneIdsNow = c(rhyGenes$geneId[1], drLimma$geneId[1], deLimma$geneId[1])
geneSymbolsNow = unname(getSYMBOL(geneIdsNow, 'org.Mm.eg.db'))

df = as_tibble(t(emat[geneIdsNow,]))
colnames(df) = geneSymbolsNow
df$sample = colnames(emat[geneIdsNow,])
df = df %>%
  inner_join(sm, by = 'sample') %>%
  select(sample, cond, time, geneSymbolsNow) %>%
  gather(key = 'geneSymbol', value = 'expr', geneSymbolsNow) %>%
  mutate(geneSymbol = factor(geneSymbol, levels = geneSymbolsNow))

ggplot(df) +
  facet_grid(geneSymbol ~ cond, scales = 'free_y') +
  geom_point(aes(x = time, y = expr, shape = cond, color = geneSymbol), size = 2) +
  labs(x = 'Zeitgeber time (h)', y = 'Expression (norm.)') +
  scale_shape(solid = FALSE, guide = FALSE) +
  scale_color_brewer(type = 'qual', palette = 'Dark2', guide = FALSE) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6))

