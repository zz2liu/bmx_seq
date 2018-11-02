input = snakemake@input
output = snakemake@output
#wildcards = snakemake@wildcards
#config = snakemake@config
#print(input)
#save.space()

master_contrast_start = 7

mastersheets = input[names(input)!='']
master = read.csv(mastersheets[[1]], row.names=1, check.names=F)
for (a in names(mastersheets)[-1]) {
    curr = read.csv(mastersheets[[a]], row.names=1, check.names=F)
    curr = curr[,master_contrast_start:ncol(curr)]
    colnames(curr) = paste0(colnames(curr), '_', a)
    master = cbind(master, curr[rownames(master),])
}

#save.image()
write.csv(master, output[[1]])
