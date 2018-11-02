#!/usr/bin/env python3
import numpy as np, pandas as pd, seaborn as sns
from matplotlib import pyplot as plt
from pathlib import Path

def parseTrimming(infile):
    """return the totalReads"""
    for line in infile:
        if line.startswith('Total reads processed:'):
            res = line.split(':')[1].strip().replace(',','')
            return int(res)

def test_parseTrimming():
    input = '''=== Summary ===

Total reads processed:              35,618,681
Reads with adapters:                11,413,864 (32.0%)
Reads written (passing filters):    35,618,681 (100.0%)'''.splitlines()
    res = parseTrimming(input)
    assert res == 35618681

def parseMapping(infile):
    """return the totalReads,uniquelyMapped, multiMapped, toomanyMapped"""
    #raw = {k.strip(): v.strip()
    #    for line in infile if ' |' in line
    #        for k, v in [line.split(' |', 1)]}
    raw = dict((x.strip() for x in line.split(' |', 1))
        for line in infile if ' |' in line)
    #print(raw)
    res = tuple(int(raw[k]) for k in (
        'Number of input reads', 'Uniquely mapped reads number', 
        'Number of reads mapped to multiple loci','Number of reads mapped to too many loci')) 
    return res

def test_parseMapping():
    input = '''
       Mapping speed, Million of reads per hour |       597.02

                          Number of input reads |       24378270
                      Average input read length |       74
                                    UNIQUE READS:
                   Uniquely mapped reads number |       17921860
                        Uniquely mapped reads % |       73.52%
        Number of reads mapped to multiple loci |       6126647
             % of reads mapped to multiple loci |       25.13%
        Number of reads mapped to too many loci |       204806
             % of reads mapped to too many loci |       0.84%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.01%
                 % of reads unmapped: too short |       0.30%
                     % of reads unmapped: other |       0.20%
'''.splitlines()
    res = parseMapping(input)
    assert tuple(res) == (24378270,17921860, 6126647,204806)

def collectData(trimLogs, mapLogs, samples):
    """return a dataframe of [sample: trimDropped, mapDropped, mapped]"""
    total = np.array([parseTrimming(open(f)) for f in trimLogs])
    trimmed, unique, multi, toomany = np.array([parseMapping(open(f)) for f in mapLogs]).T
    #import pdb; pdb.set_trace()
    res = pd.DataFrame(np.array([total-trimmed, trimmed-multi-toomany-unique, multi+toomany, unique]).T, 
        columns=['trimmDropped','notMapped','multiMapped','uniquelyMapped'], index=samples)
    return res

def test_collectData():
    samples = [f'Sample_EW2{i}' for i in [4,5,6]]
    trimLogs = [f'tests/data/{sample}.R1.fastq.gz_trimming_report.txt' for sample in samples]
    mapLogs = [f'tests/data/{sample}.Log.final.out' for sample in samples]
    res = collectData(trimLogs, mapLogs, samples)
    #print (res)
    assert res.shape == (3,4)

def _plotData(data, ylabel='', outfile=None):
    sns.set()
    ax = data.plot.bar(stacked=True)
    ax.set(xlabel='Sample', ylabel=ylabel)
    plt.legend(bbox_to_anchor=(1,1))
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, bbox_inches='tight')

def plotData(data, outfile):
    #import pdb; pdb.set_trace()
    data = data.iloc[:, ::-1]
    _plotData(data, 'Reads', outfile)

def plotDataPerc(data, outfile=None):
    data = data.iloc[:, ::-1]
    perc = (data.T/data.sum(1)).T
    _plotData(perc, 'Percentage of Reads', outfile)

def main(snakemake):
    data = collectData(snakemake.input.trimming, snakemake.input.mapping, snakemake.params.samples)
    data.to_csv(snakemake.output.csv)
    plotData(data, snakemake.output.svg)
    plotDataPerc(data, snakemake.output.percSvg)

def test_main():
    class obj:
        class input:
            trimming=[f'tests/data/Sample_EW2{i}.R1.fastq.gz_trimming_report.txt' for i in [4,5,6]]
            mapping=[f'tests/data/Sample_EW2{i}.Log.final.out' for i in [4,5,6]]
        class params:
            samples=[f'Sample_EW2{i}' for i in [4,5,6]]
        class output:
            csv='/tmp/mappingSummary.csv'
            svg='/tmp/mappingSummary.svg'
            percSvg='/tmp/mappingSummaryPerc.svg'
    main(obj)

if __name__ == '__main__':
    main(snakemake)
