#take args
input = snakemake@input #count,info
output = snakemake@output #csv,pdf
#print(input)
#print(output)
# read and parse input
count = read.csv(input$count, row.names=1, check.names=F)
info = read.csv(input$info, row.names=1, check.names=F)
save.image()

tmp = aggregate(count, by=list(geneType=info$gene_type), FUN=sum)
rownames(tmp) = tmp[,1]; tmp[,1]=NULL
#gene.summary = tmp

# sort by max percentage in all samples, to catch the outlier sample
prop = t(t(tmp)/colSums(tmp))
ord = order(apply(prop, 1, max), decreasing=T)
gene.summary = tmp[ord,]

# write csv
write.csv(gene.summary, output$csv)

# write pdf
tmp = as.matrix(gene.summary)
ROW_MAX=7
# merge genetypes with lower percentages to others
if (nrow(tmp)>ROW_MAX+1) {
    tmp = rbind(tmp[1:7,], other=colSums(tmp[7:nrow(tmp),]))
}

doPlot = function(x, outpdf) {
    par(mar=c(6,4,2,2))
    pdf(outpdf, width=12)
    barplot(x, col=1:nrow(x) + 1,
        las=2, legend=T, xlim=c(0,ncol(x)*2),
        ylab='Total reads of GeneTypes')
    prop = t(t(x)/colSums(x))
    barplot(prop, col=1:nrow(prop) + 1,
        las=2, legend=T, xlim=c(0,ncol(prop)*2),
        ylab='Proportion of GeneTypes')
    dev.off()
}

doPlot(tmp, output$pdf)
