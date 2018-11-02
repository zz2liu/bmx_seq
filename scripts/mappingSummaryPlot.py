#!/usr/bin/env python3
## not working because of the young stage of ggplot
from rpy2gt
def ggplotData(data, outfile):
    df = melt(as.matrix(data))
    colnames(df)=c('Sample','Cat', 'Reads')
    save.image()
    ggplot(df, aes(fill=Cat, y=Reads, x=Sample)) +
        geom_bar( stat="identity") + # coord_flip()
        scale_y_continuous(expand=c(0,0)) + #labels = function(x) paste0(x*100, "%")
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5))

def test_plotData():
    data = readData('tests/data/mappingSummary.csv')
    tmp = pd.melt(data.reset_index(), ['index'], var_name='Cat', value_name='Reads')
    df = tmp.rename(columns=dict(index='Sample', inplace=True))
    g = ggplot(df, aes('Sample', 'Reads', fill='Cat'))
        geom_bar(stat='identity')
    
def readData(infile):
    res = pd.read_csv(infile, index_col=0)
    return res

def main(snk):
    data = readData(snk.input.csv)
    plotData(data, snk.output.pdf)

if __name__ == '__main__':
    main(snakemake)
