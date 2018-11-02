### read and align geneCount(a tsv) and geneInfo
gene_norm_rpkm = function(geneCountFile, geneInfoFile) {
    # read geneCount,geneInfo
    geneCount = read.csv(geneCountFile, row.names=1, check.names=F)
    geneInfo = read.csv(geneInfoFile, row.names=1, check.names=F)
    stopifnot(rownames(geneCount)==rownames(geneInfo))
    
    # extract geneK from geneInfo
    geneK = geneInfo$medianTxLength / 1000

    # calc sampleM from geneCount
    sampleM = colSums(geneCount) / 1000000

    # calc rpkm
    res = t(t(geneCount)/sampleM) / geneK
    return(res)
}

## input$info,count > output$info,count
input = snakemake@input
output = snakemake@output
geneCountFile = input$count #'hgmmStar.hg/geneCount.txt'
geneInfoFile = input$info #'~/local/data/hg38/geneInfo.csv'

save.image()
rpkm = gene_norm_rpkm(geneCountFile, geneInfoFile)

write.csv(rpkm, output[[1]])

