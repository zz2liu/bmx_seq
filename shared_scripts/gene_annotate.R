### read and align geneCount(a tsv) and geneInfo
gene_annotate = function(geneCountFile, geneInfoFile) {
    #read geneCount and filter genes
    raw = read.delim(geneCountFile, row.names=1, check.names=F)
    sele = rowSums(raw)>0
    geneCount = raw[sele,]

    #read geneInfo,align with geneCount,name rows with name|id.version
    raw = read.csv(geneInfoFile, row.names=1, check.names=F)
    geneInfo = raw[row.names(geneCount),]
    row.names = paste(geneInfo$gene_name, rownames(geneInfo), sep='|')
    rownames(geneInfo) = row.names
    rownames(geneCount) = row.names

    return (list( count=geneCount, info=geneInfo ))
}

.writeTables = function(tables, outDir) {
    #write fixed tables
    dir.create(outDir, recursive=T)
    for (i in names(tables)) {
        write.csv(tables[[i]], paste0(outDir, '/', i, '.csv'))
    }
}

## input$info,count > output$info,count
input = snakemake@input
output = snakemake@output
geneCountFile = input$count #'hgmmStar.hg/geneCount.txt'
geneInfoFile = input$info #'~/local/data/hg38/geneInfo.csv'
save.image()

fixed = gene_annotate(geneCountFile, geneInfoFile)
write.csv(fixed$info, output$info)
write.csv(fixed$count, output$count)



