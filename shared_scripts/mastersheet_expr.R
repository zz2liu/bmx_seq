input = snakemake@input
output = snakemake@output

# read the input files.
mastersheet = read.csv(input$mastersheet, row.names=1, check.names=F)
rpkm = read.csv(input$rpkm, row.names=1, check.names=F)
vst = read.csv(input$vst, row.names=1, check.names=F)

# write the aligned expr files
write.csv(rpkm[rownames(mastersheet),], output$rpkm)
write.csv(vst[rownames(mastersheet),], output$vst)


