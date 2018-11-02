input = snakemake@input
output = snakemake@output

# read and parse input
data = read.csv(input$csv, row.names=1, check.names=F)

# write pdf
#x = t(as.matrix(data))
#x = x[nrow(x):1,]
pdf(output$pdf, width=12)
#par(mar=c(12,6,2,2))
#barplot(x, col=1:nrow(x) + 1,
#    las=2, legend=T, xlim=c(0,ncol(x)*2),
#    ylab='Total Reads')
#prop = t(t(x)/colSums(x))
#barplot(prop, col=1:nrow(prop) + 1,
#    las=2, legend=T, xlim=c(0,ncol(prop)*2),
#    ylab='Proportion of Reads')

# ggplot
require(reshape2)
require(ggplot2)
require(scales)
df = melt(as.matrix(data))
colnames(df)=c('Sample','Cat', 'Reads')
#save.image()
ggplot(df, aes(fill=Cat, y=Reads, x=Sample)) +
    geom_bar( stat="identity") + # coord_flip()
    scale_y_continuous(expand=c(0,0)) + #labels = function(x) paste0(x*100, "%")
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5))
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))

# Stacked Percent
ggplot(df, aes(fill=Cat, y=Reads, x=Sample)) +
    geom_bar( stat="identity", position="fill") +
    scale_y_continuous(labels = scales::percent, expand=c(0,0)) + #labels = function(x) paste0(x*100, "%")
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5))
dev.off()
