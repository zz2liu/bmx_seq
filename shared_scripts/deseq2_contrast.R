input = snakemake@input
output = snakemake@output
params = snakemake@params #padj_threshold, fc_threshold
wildcards = snakemake@wildcards
groupCol = 'group'


# parse contrast
cc = strsplit(wildcards$contrast, split='-')[[1]]
i = wildcards$contrast
print (i)

# load dds and calc diffexpr
suppressMessages(suppressWarnings(
    require('DESeq2')
))
dds = readRDS(input$rds) #DESeq(dds) #todo: fitType warning to be documented
res = results(dds, contrast=c(groupCol, cc))

# plotMA
pdf(output$maplot) #paste0(i, '.plotMA.pdf'))
DESeq2::plotMA(res, alpha=.05, ylim=c(-3,3), main=i)
#identify(res$baseMean, res$log2FoldChange, labels=rownames(res))
dev.off()

# parse group and expr
group = read.csv(input$sampleInfo, row.names=1, check.names=F)[,groupCol]
expr = read.csv(input$vst, row.names=1, check.names=F)

# add expr to DE res
x = as.data.frame(res)
x = x[order(x$stat),]
xx = cbind(x, expr[rownames(x), c(which(group==cc[1]), which(group==cc[2]))])

## write DE files
curr = as.data.frame(xx)
#independent filtering passed and failed
mask = is.na(curr$padj)
write.csv(curr[!mask,], output$accepted) #paste0(i, '.StatWithVst.accepted.csv'))
write.csv(curr[mask,], output$rejected) #paste0(i, '.StatWithVst.rejected.csv'))
#significant by padj
sele = which(curr$padj< params$padj_threshold & abs(curr$log2FoldChange) >= log2(params$fc_threshold))
# write the csv even if there is nothing significant, to satisfy the output
if (length(sele)) {
    write.csv(curr[sele,], output$sig)
    with(curr[sele,], print (table(sign(stat))))
} else {
    file.create(output$sig)
}
    

## write rnk files for GSEA
#writeRnk(output$accepted)
