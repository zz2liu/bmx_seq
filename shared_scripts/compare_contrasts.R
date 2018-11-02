input = snakemake@input
output = snakemake@output
params = snakemake@params
if(is.null(params$contrastComparisons)) {
    params$contrastComparisons = params$contrasts
}

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

## summarize plot using SuperExactTest
maxpadj = params$padj_threshold
minfc = params$FC_threshold
total = 20000
genesets = list(up=list(), down=list(), both=list())
for (i in names(all_contrasts)) {
    curr = as.data.frame(all_contrasts[[i]])
    up = with(curr, which(padj<=maxpadj & log2FoldChange>=log2(minfc)))
    down = with(curr, which(padj<=maxpadj & log2FoldChange<=-log2(minfc)))
    both = with(curr, which(padj<=maxpadj & abs(log2FoldChange)>=log2(minfc)))
    genes = sub('(.+)\\|.+', '\\1', rownames(curr), perl=T)
    genesets$up[[i]] = unique(genes[up])
    genesets$down[[i]] = unique(genes[down])
    genesets$both[[i]] = unique(genes[both])
}

if (length(genesets$both) == 1) { #only one contrast
    for (f in output) {
        system(paste0('touch ', f))
    }
    quit()
}

#require('SuperExactTest')
#save.image()
#for (i in names(genesets)) {
#    data = genesets[[i]]
#    # shortcut for only one contrast with any significant genes
#    sele = sapply(data, function(x) {
#        length(x)>0 })
#    if (sum(sele) < 2) {
#        system(paste0('touch ', output[[i]]))
#        next
#    }
#
#    res = supertest(data, n=20000)
#    write.csv(summary(res)$Table, file=output[[i]])
#=======
indicateSets = function(sets) {
    row.names = Reduce(union, sets)
    res = matrix(0, nrow=length(row.names), ncol=length(sets))
    rownames(res) = row.names
    colnames(res) = names(sets)
    for (i in names(sets)) { res[sets[[i]], i]=1}
    res
}
#tmp = indicateSets(data)
require(UpSetR)
require(limma)
save.image()
for (i in names(genesets)) {
    ind = as.data.frame(indicateSets(genesets[[i]]))
    write.csv(ind, file=output[[i]])
    pdf(output[[paste0(i, '_pdf')]], width=20, onefile=F)
    for (cmp in params$contrastComparisons) {
        dat = ind[cmp]
        if (ncol(dat) <= 3) {
            vennDiagram(dat)
        } else {
            upset(dat, ncol(dat), NA, order.by='freq', decreasing=T)
        }
    }
    dev.off()
}

############################################################## old_par
#save.image()
#quit()
#require('SuperExactTest')
#for (i in names(genesets)) {
#    res = supertest(genesets[[i]], n=20000)
#    write.csv(summary(res)$Table, file=output[[i]])
#
#    pdf(output[[paste0(i, '_pdf')]], width=15)
#    old_par = par()
#    #par(mfrow=c(2,1), mar=c(5,8,4,1))
#    par(mar=c(5,15,4,1))
#    #plot(res, sort.by="size", degree=2:length(res)) #, main=paste(i, 'regulated genes'))
#    plot(res, Layout='landscape', degree=2:length(res)) #, main=paste(i, 'regulated genes'))
#    dev.off()
#    par(old_par)
#}

