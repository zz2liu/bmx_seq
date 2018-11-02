input = snakemake@input
params = snakemake@params
output = snakemake@output
#SAMPLE_INFO = snakemake@config$SAMPLE_INFO
#
#sampleInfo = read.csv(SAMPLE_INFO, row.names=1)
#sampleInfo$sample = rownames(sampleInfo)

# start
#require(DESeq2)
require(ggplot2)
require(ggrepel)

#sizevst = readRDS(input$vst)
vst = read.csv(input$vst, row.names=1, check.names=F)
sampleInfo = read.csv(input$sampleInfo, row.names=1, check.names=F)
sampleInfo$sample = rownames(sampleInfo)

#save.image()
tmp = prcomp(t(vst))
x = tmp$x
.var = tmp$sdev ^ 2
percentVar = round(100 * .var/sum(.var))
con = file(output$csv, open='wt')
writeLines(paste(c('#percentVar', percentVar), collapse=','), con)
write.csv(x, con)
close(con)

data = cbind(x, sampleInfo)
pdf(output$pdf) #'plotPca_top500.pdf')
p = ggplot(data, aes_string('PC1', 'PC2', 
    shape=params$shape, color=params$color, label=params$label, size=params$size))
p = p + geom_point(size=params$pointSize, alpha=params$pointAlpha) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))# +
if (params$equalScale) {
    p = p + coord_fixed()
}
if (!is.null(params$label)) {
    p = p + geom_text_repel(size=params$labelSize, color=params$labelColor)
}
plot(p)
dev.off()


