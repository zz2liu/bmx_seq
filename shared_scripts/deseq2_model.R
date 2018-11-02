# DESeq2 preprocess, normalization, modeling
# @input
input = snakemake@input
# @params
params = snakemake@params
# @output
output = snakemake@output

# funs
assert = function(test, msg) {
    if (!all(test)) stop(msg)
}
# parse contrasts list from the commandline, check againt sampleInfo later
contrasts = lapply(params$contrasts, function(x) {strsplit(x, split='-')[[1]]})
names(contrasts) = params$contrasts

# read sampleInfo
# todo: move this to config
groupCol = params$group; batchCol=params$batch
#raw = read.csv(input$sampleInfo, row.names=1, header=F, skip=1, as.is=T)
sampleInfo = read.csv(input$sampleInfo, row.names=1, as.is=T)
    # as.is is used to avoid factor trimming problems
##############################################################
#save.image()
assert(groupCol %in% colnames(sampleInfo),
    'InputError: the groupCol not found in the sampleInfo')
## check contrast units
assert(unlist(contrasts) %in% sampleInfo[,groupCol],
    paste('contrasts: some sampleGroups not found in your sampleInfo:',
        paste(setdiff(unlist(contrasts), sampleInfo[,groupCol]), collapse=',')))
## filter sampleInfo by contrast units
sele = sampleInfo[,groupCol] %in% unlist(contrasts)
sampleInfo = sampleInfo[sele,,drop=F] #drop necessary, otherwise get a vector
sampleInfo[,1] = as.factor(sampleInfo[,1])

if (!is.null(batchCol)) {
    assert(batchCol %in% colnames(sampleInfo),
        'InputError: the batchCol not found in the sampleInfo')
    sampleInfo[,2] = as.factor(sampleInfo[,2])
}

# read geneCoutn, Align and filter to sampleInfo, keep only nonzero genes
#todo: readGeneCount(countFile, sampleInfo=NULL)
geneCount = read.csv(input$geneCount, 
    row.names=1, check.names=F)
## check sampleInfo
assert(rownames(sampleInfo) %in% colnames(geneCount),
    paste('sampleInfo: some sampleNames not found in geneCount:', 
        paste(setdiff(rownames(sampleInfo), colnames(geneCount)), collapse=',')))
## align to sampleInfo
countNz = geneCount[rowSums(geneCount)>0, rownames(sampleInfo)]

## write the filtered input
write.csv(sampleInfo, output$sampleInfo) #'sampleInfo.filtered.csv')
write.csv(countNz, output$geneCount) #'geneCount.filtered.csv')

## prepare dataset
if (!is.null(batchCol)) {
    design = as.formula(paste0('~', batchCol, '+', groupCol))
} else {
    design = as.formula(paste0('~', groupCol))
}
suppressMessages(suppressWarnings(
    require('DESeq2')))
dds = DESeqDataSetFromMatrix(countData=countNz, colData=sampleInfo, design=design)

## vst_blind and pca plot
#todo: deseqSummary(dds)
vst = varianceStabilizingTransformation(dds, blind=T) #fitType is quickfix to suppress warnings
saveRDS(vst, output$vstRds)
expr = assay(vst)
write.csv(expr, output$vst) #'expr.vst_blind.csv')


## fit the model 
dds = DESeq(dds) #todo: fitType warning to be documented
saveRDS(dds, output$rds) #'dds.rds')


