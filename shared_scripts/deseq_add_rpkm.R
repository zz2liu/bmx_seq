## add rpkm to de csv
deseq_add_rpkm = function(deFile, rpkmFile) {
    # read and check
    rpkm = read.csv(rpkmFile, row.names=1, check.names=F)
    de = read.csv(deFile, row.names=1, check.names=F)
    sampleNames = colnames(de)[7:ncol(de)]
    stopifnot(sampleNames %in% colnames(rpkm))
    stopifnot(rownames(de) %in% rownames(rpkm))
    
    # extract rpkm and collate
    tmp = rpkm[rownames(de), sampleNames]
    log2Rpkm = log2(tmp+1)
    value = rowMeans(log2Rpkm)

    res = cbind(avgLog2Rpkm=value, de)
    return(res)
}

deseq_plot_rpkm = function(deRpkm)
{
    tmp = as.data.frame(data.matrix(deRpkm[, c('avgLog2Rpkm', 'log2FoldChange')]))
    with(tmp,
        plot(log2FoldChange ~ avgLog2Rpkm)
    )
}

## input$info,count > output$info,count
input = snakemake@input
output = snakemake@output
deFile = input$diffexpr
rpkmFile = input$rpkm

# create empty output if input is empty
if (file.info(deFile)$size==0) {
    file.create(output$csv)
    file.create(output$png)
} else {
    # calc and save res
    res = deseq_add_rpkm(deFile, rpkmFile)
    write.csv(res, output[['csv']])

    # plot and export to png
    png(output[['png']])
    deseq_plot_rpkm(res)
    dev.off()
}


