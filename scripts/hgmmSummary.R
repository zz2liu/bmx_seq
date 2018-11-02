#take args
input = snakemake@input #mm,hg
output = snakemake@output #csv,pdf

# read and parse input
mm.summary = read.delim(input$mm, row.names=1, check.names=F)
hg.summary = read.delim(input$hg, row.names=1, check.names=F)
#save.image()

rownames(mm.summary) = paste0('mm.', rownames(mm.summary))
rownames(hg.summary) = paste0('hg.', rownames(hg.summary))
hgmm.summary = rbind(hg.summary, mm.summary)
#tmp = sub('Assigned', 'oneGene', sub('Unassigned_NoFeatures','noGene', sub('Unassigned_Ambiguity', 'manyGenes', rownames(hgmm.summary))))
#rownames(hgmm.summary) = tmp
hgmm.summary=hgmm.summary[rowSums(hgmm.summary)>0,]

# write csv
write.csv(hgmm.summary, output$csv)

# write pdf
x = as.matrix(hgmm.summary)
pdf(output$pdf, width=12)
par(mar=c(12,4,2,2))
barplot(x, col=1:nrow(x) + 1,
    las=2, legend=T, xlim=c(0,ncol(x)*2),
    main='Total Reads')
prop = t(t(x)/colSums(x))
barplot(prop, col=1:nrow(prop) + 1,
    las=2, legend=T, xlim=c(0,ncol(prop)*2),
    ylab='Proportion of Reads')
dev.off()


