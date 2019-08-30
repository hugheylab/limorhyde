library('foreach')


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
  mappingProbeIntermediate = featureDf[!is.na(featureDf[[interName]]) & featureDf[[interName]] != '',
                                       c('ID', interName)]
  colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')

  mapTmp1 = eval(parse(text = sprintf('%s.db::%s', substr(dbName, 1, 9), dbName)))
  mapTmp2 = AnnotationDbi::mappedkeys(mapTmp1)
  mapTmp3 = as.list(mapTmp1[mapTmp2])
  geneId = do.call(c, mapTmp3)

  geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3),
                                sapply(mapTmp3, length), SIMPLIFY = FALSE))
  if (dbName == 'org.Hs.egUNIGENE2EG') {
    geneInter = sub('Hs.', '', geneInter, fixed = TRUE)}
  mappingIdInter = data.frame(geneId, geneInter, stringsAsFactors = FALSE)
  return(merge(mappingIdInter, mappingProbeIntermediate, by = 'geneInter', sort = FALSE))}


calcExprsByGene = function(eset, mapping) {
  geneIds = unique(mapping$geneId)
  exprsByGene = matrix(nrow = length(geneIds), ncol = ncol(eset),
                       dimnames = list(geneIds, Biobase::sampleNames(eset)))
  for (geneId in geneIds) {
    exprsTmp = Biobase::exprs(eset)[mapping$probeSet[mapping$geneId == geneId],, drop = FALSE]
    if (nrow(exprsTmp) == 1) {
      exprsByGene[geneId,] = exprsTmp
    } else {
      exprsByGene[geneId,] = matrixStats::rowMedians(t(exprsTmp), na.rm = TRUE)}}
  return(exprsByGene)}


fixGeoSampleNames = function(sampleNames) {
  sampleNames = paste0(toupper(sampleNames), '_')
  regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
  sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength - 1),
                          sampleNames, attr(regexResult, 'match.length'))
  return(unname(sampleNamesNew))}


getGeneMapping = function(eset){
  featureDf = pData(eset)
  idx = sapply(featureDf, is.factor)
  featureDf[idx] = lapply(featureDf[idx], as.character)
  featureDf$RefSeq = sapply(featureDf$GB_ACC, function(x) strsplit(x, split = '.', fixed = TRUE)[[1]][1])
  return(getGeneProbeMappingAnno(featureDf, dbName = 'org.Mm.egREFSEQ2EG', interName = 'RefSeq'))}
