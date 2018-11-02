# DESeq2 likelihood ratio test (LRT)
# @input
input = snakemake@input
# @output
output = snakemake@output

# read sampleInfo
dds = readRDS(input$rds)

## prepare dataset
suppressMessages(suppressWarnings(
    require('DESeq2')))

## later: reduced = paste0(format(design(dds)), ' - group')
if (format(design(dds))=="~group") { #only model with group
    ddsLrt = DESeq(dds, test='LRT', reduced =~ 1)
} else {
    ddsLrt = DESeq(dds, test='LRT', reduced =~ batch)
}
res = results(ddsLrt)
ord = order(-res$stat)
#save.image() # debug
write.csv(res[ord, c('stat', 'pvalue', 'padj')], output$csv)

