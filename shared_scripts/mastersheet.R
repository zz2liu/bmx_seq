input = snakemake@input
output = snakemake@output
params = snakemake@params

suppressMessages(suppressWarnings(
    require('DESeq2')
))
dds = readRDS(input$rds)

# collect contrast results
all_contrasts = list()
sele = c('log2FoldChange', 'pvalue', 'padj')
for (i in params$contrasts) {
    cc = strsplit(i, split='-')[[1]]
    res = results(dds, contrast=c('group', cc))
    all_contrasts[[i]] = res[,sele]
}
contrasts_res = Reduce(cbind, all_contrasts)
colnames(contrasts_res) = paste0(
    rep(names(all_contrasts), each=length(sele)),
    '.',
    rep(sele, times=length(all_contrasts)))

# collect anova results
reduced = paste0(format(design(dds)), ' - group')
ddsLrt = DESeq(dds, test='LRT', reduced = as.formula(reduced))
res = results(ddsLrt)
sele = c('stat', 'pvalue', 'padj')
anova_res = res[,sele]
colnames(anova_res) = paste0('anova', '.', sele)
# order of the genes in mastersheet
ord = order(-res$stat)
genes = rownames(res)[ord]

# average RPKM of all the samples included
sampleInfo = read.csv(input$sampleInfo, row.names=1, check.names=F)
tmp = read.csv(input$rpkm, row.names=1, check.names=F)
rpkm = tmp[,rownames(sampleInfo)]
avgLog2Rpkm = rowMeans(log2(rpkm+1))

# construst the mastersheet
gene_name = sub('(.+)\\|.+', '\\1', genes, perl=T)
gene_id = sub('.+\\|([^.]+)\\..+', '\\1', genes, perl=T)
mastersheet = data.frame(gene_name=gene_name, gene_id=gene_id,
    avgLog2Rpkm=avgLog2Rpkm[genes], 
    anova_res[genes,], 
    contrasts_res[genes,])
write.csv(mastersheet, output$csv)

