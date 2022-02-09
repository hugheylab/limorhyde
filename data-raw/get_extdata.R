library('data.table')

parentDir = 'data-raw'
outputDir = file.path('inst', 'extdata')
study = 'GSE34018'

seeker::seekerArray(list(study = study, geneIdType = 'entrez'), parentDir)

metadata = fread(file.path(parentDir, study, 'sample_metadata.csv'))
metadata = metadata[, .(sample_id, title, genotype = `genotype/variation:ch1`)]
metadata[, time := as.numeric(sub('_.*', '', sub('.*_ZT', '', title)))]
metadata[, cond := factor(
  genotype, c('wild-type', 'Reverb alpha/beta DKO'), c('wild-type', 'knockout'))]
setorderv(metadata, c('cond', 'time'))

emat = qs::qread(file.path(parentDir, study, 'gene_expression_matrix.qs'))
emat = emat[, metadata$sample_id]
emat = log2(emat)

qs::qsave(metadata, file.path(outputDir, paste0(study, '_metadata.qs')))
qs::qsave(emat, file.path(outputDir, paste0(study, '_expression_data.qs')))

unlink(file.path(parentDir, study), recursive = TRUE)
