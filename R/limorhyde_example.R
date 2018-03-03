# create an example data frame
nSamples = 12
df = data.frame(
   sample = paste0('sample_', 1:nSamples),
   genotype = factor(rep(c('WT', 'KO'), each=nSamples/2), levels=c('WT', 'KO')),
   zt = rep(seq(0, 24-24/nSamples*2, 24/nSamples*2), times=2),
   stringsAsFactors=FALSE)

# call limorhyde and examine its output
limoResult = limorhyde(df, 'zt')
limoResult$df
limoResult$timeColnames

# create a design matrix that could be used with methods such as limma,
# in order to detect differentially rhythmic genes
design = model.matrix(~ genotype * (zt_cos + zt_sin), data=limoResult$df)
