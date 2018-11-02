"""
USAGE:
export PATH=/home/zl99/anaconda3/bin:/home/zl99/code/ngs/pipelines:$PATH
alias snakeSlurm="snakemake --cores 99 --cluster-config cluster.json \
    --cluster 'sbatch --mem {cluster.mem} -c {threads} -J {rule} -o {rule}.%j.out'"
snakeSlurm -pn
snakeSlurm -p
alias snakeDag="snakemake --forceall --dag | dot -Tpng > dag.png"

# References:
- snakemake deployment: http://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
- snakemake tutorial: http://slowkow.com/notes/snakemake-tutorial/
- bash strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
"""
# later: add software requirements --use-conda
# later: use wrapper for star alignment
# later: report the hgmmStar, (bamstat, genetype on hg/mm), (pca, degs on hg/mm)
# later: hgmmSummary include nRaw, nTrimmed
# later: star: paired-ended, two passes

import sys, textwrap, pandas as pd, snakemake as smk
import collections
from os import path

configfile: 'config.yaml'
#config['SAMPLE_INFO'] = config['SAMPLE_FILE'] #for back compatibility
# Globals
# ============================================================  
DATA_DIR = config['DATA_DIR']
SNAKE_DIR = config['SNAKE_DIR']
GENOMES = config['GENOMES']
INPUT_DIR = config['INPUT_DIR']
SAMPLE_FILE = config.get('SAMPLE_FILE', config['SAMPLE_INFO']) #quick fix

# parse samples for sampleFile, or make a sampleFile and exit
if not path.isfile(SAMPLE_FILE): 
    smk.shell("ls -d {INPUT_DIR}/* | sort -V " #sort -V in version mode.
    "| sed -E 's%{INPUT_DIR}/(.*)%\1,A,a%' "
    "| sed '1isample,group,batch' > {SAMPLE_FILE}")
    raise SystemExit('sample file not found: a new file is created based on the input directory structure, check and edit it before run this workflow again')

# import pdb; pdb.set_trace()
df = pd.read_csv(SAMPLE_FILE, index_col=0)
SAMPLES = [str(x) for x in df.index]

ALL_TRIMMED = expand("trim_galore/{sample}.R1_trimmed.fq.gz", sample=SAMPLES)
ALL_MERGED = expand("catFastq/{sample}.R1.fastq.gz", sample=SAMPLES)
all_rpkms = [f"{genome}/sorted/{sample}.rpkm.bw" for genome in GENOMES for sample in SAMPLES]
_ALL_ANALYSES = [f"{genome}.deseq2.{analysis}"
    for genome in GENOMES #config['deseq2'].keys()
        for analysis in config['deseq2'].get(genome, [])
    ]

_ALL_CONTRASTS = [f"{genome}.deseq2.{analysis}/{contrast}"
    for genome in GENOMES #config['deseq2'].keys()
        for analysis in config['deseq2'].get(genome, [])
            for contrast in config['deseq2'][genome][analysis]['contrasts']
    ]
#import pdb; pdb.set_trace()

# rules
# ============================================================ 
#rule test:
#    input: all_rpkms

rule report:
    input:
        'hgmmSummary/hgmmSummary.pdf',
        expand("{genome}.gene_count/geneSummary.{ext}", genome=config['GENOMES'], ext=['csv', 'pdf']),
        [f"{prefix}.StatWithVst.padj_05.addRpkm.csv" for prefix in _ALL_CONTRASTS],
        [f"{anadir}/mastersheet.rpkm.csv" for anadir in _ALL_ANALYSES],
        [f"{anadir}/compare_contrasts.up.csv" for anadir in _ALL_ANALYSES],
        [f"{genome}.deseq2/mastersheet.csv" for genome in config['GENOMES']],
        pca=[f"{anadir}/plotPca.pdf" for anadir in _ALL_ANALYSES],
        #all_rpkms
    output:
        "report.html"
    script:
        "../shared/scripts/report.py"

    #input: all_rpkms
TEST_SAMPLES = 'Sample_EW16 Sample_EW17 Sample_EW18 Sample_EW19 Sample_EW20 Sample_EW21 Sample_EW1 Sample_EW2 Sample_EW3'.split()
rule test:
    input:
        hg_rpkms = [f"{genome}/sorted/{sample}.rpkm.bw" for genome in ['hg38'] for sample in TEST_SAMPLES]
        #[f"{prefix}.StatWithVst.padj_05.addRpkm.csv" for prefix in _ALL_CONTRASTS],

rule all_deseqModel:
    input:
        pca=[f"{anadir}/plotPca.pdf" for anadir in _ALL_ANALYSES],

rule all_deseqDeg:
    input:
        degs=[f"{prefix}.StatWithVst.padj_05.csv" for prefix in _ALL_CONTRASTS],


include: "../shared/rules/input.rules" #cat and trim
include: "rules/alignment.rules" #star align to hgmm, extract_hg and extract_mm
include: "../shared/rules/counts.rules" # count, annotate, summary, norm
include: "../shared/rules/diffexpr.rules" #deseq2: model, contrast, compare, mastersheet
include: "../shared/rules/deepTools.rules" #bamcoverage
include: "../shared/rules/samtools.rules" #index
#include: path.join(SNAKE_DIR, "rules/alignment.rules") #align to hgmm, extract_hg and extract_mm
    #unecessary to use SNAKE_DIR, if point the snakefile here.

def _find_analyses_under_genome(wildcards):
    return collections.OrderedDict(
        (f'{ana}',   wildcards.genome + f'.deseq2.{ana}/mastersheet.csv')
        for ana in config['deseq2'][wildcards.genome])

#rule checkin_config:
#    output:
#    script:
rule merge_mastersheets:
    input:
        unpack(_find_analyses_under_genome)
    output:
        '{genome}.deseq2/mastersheet.csv'
    script:
        '../shared/scripts/merge_mastersheets.R'

# add a column of mean RPKM of all samples in a deseq2 output
rule deseq_add_rpkm:
    input:
        diffexpr="{genome}.deseq2.{analysis}/{contrast}.StatWithVst.padj_05.csv",
        rpkm=rules.gene_norm_rpkm.output[0]
    output:
        csv="{genome}.deseq2.{analysis}/{contrast}.StatWithVst.padj_05.addRpkm.csv",
        png="{genome}.deseq2.{analysis}/{contrast}.StatWithVst.padj_05.addRpkm.png"
    script:
        "../shared/scripts/deseq_add_rpkm.R"

# {genome} to be implemented later using input: unpack(lambda wildcards:...)
# unmapped and unpassed reads to be added later from hgmmStar output
rule mappingSummary:
    input:
        trimming=[f'trim_galore/{sample}.R1.fastq.gz_trimming_report.txt'
            for sample in SAMPLES],
        mapping=[f'hgmmStar/{sample}/Log.final.out'
            for sample in SAMPLES],
    output:
        csv="mappingSummary/mappingSummary.csv",
        svg="mappingSummary/mappingSummary.svg",
        percSvg="mappingSummary/mappingSummaryPerc.svg"
    params:
        samples=SAMPLES
    script:
        'scripts/mappingSummary.py'

rule mappingSummaryPlot:
    input:
        csv=rules.mappingSummary.output.csv
    output:
        pdf="mappingSummary/mappingSummary.pdf"
    script:
        'scripts/mappingSummaryPlot.R'
       
rule hgmmSummary:
    input:
        hg="hg38.gene_count/featureCounts.summary",
        mm="mm10.gene_count/featureCounts.summary"
    output:
        csv="hgmmSummary/hgmmSummary.csv",
        pdf="hgmmSummary/hgmmSummary.pdf"
    script:
        path.join(SNAKE_DIR, "scripts/hgmmSummary.R")
        #"scripts/tmp.R"
