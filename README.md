# bmx_seq
BMX-Seq: Brain Metastasis Xenograft Sequencing Workflow

This is a RNA-Seq workflow mapping reads to human and mouse genome/transcriptome.

Summray of the Workflow:
1. Raw reads were trimmed off the adaptor sequences using trim_galore.
2. Reads were mapped to combined genome and transcriptome of human and mouse using STAR.
3. The alignment to each genome were counted to annotated genes seperately using subreads::featureCounts.
4. Differential expressed genes were then detected seperately for each species and comparison using DESeq2.

A example presentation of the workflow:

{TBA}

## Authors
- Zongzhi Z Liu
- Emily Wingrove

## Installation
If you just plan to use this workflow, download and extract the latest
release. If you intend to modify and further develop this workflow, fork this
reposity. Please consider providing any generally applicable modifications via
a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits
to the authors by citing the URL of this repository and, once available the
citation to our publication.

### Requirements
- Anaconda
- Snakemake
- All the other dependents will be automatically deployed into isolated environments by Snakemake.

## Usage
### Configure
Configure the workflow via editing the following files:

- sampleInfo.csv: the sampleInfo to map sampleName to group
  and batch

- config.yaml: the configurantion file for the workflow.
  Including genome and transcriptome sources, contrasts for each DESeq2
  analysis, and parameters for other steps of the workflow.

You are suggested to copy the two config files to a working directory with your input/ fastq files, then execute the workflow with an --directory pointing to it. See an example the the demo directory.

Test your configuration by performing a dry-run via
```
snakemake --use-conda -n
```
### Execute
Execute the workflow locally using $N cores via
```
snakemake --use-conda --cores $N
```

Alternatively, excecute the workflow in cluster or cloud environments (see [the Snakemake docs](http://snakemake.readthedocs.io/en/stable/executable.html) for details). You might find the cluster.json file helpful, in this case.


### Interetation of the result files in summary (see details in the generated report.html).
- trim_galore/
  - \*.fq.gz: the raw reads with adapter and poor quality basepairs trimmed off.
  - \*\_fastqc.html: the FASTQC quality analyses.

- hgmmSummary/
  - hgmmSummary.csv|pdf: A summary and plot of the mapping workflow.

- hg38/
  - \*.bam: the alignment file for human genome.

- hg38.gene_count/

  - geneCount.csv: a gene count matrix of gene_name|gene_id.version x sampleName
  - geneRpkm.csv: a Reads Per Kbp per Million total mapped (RPKM) matrix.

- hg38.deseq2.\*/: The DESeq2 analyses of all the contrasts you defined in the config.yaml

  - plotPca_top500.pdf: the PCA plot using the top 500 genes with large variance.
  - expr.vst_blind.csv: the count matrix normalized with DESeq2
    Variance Stabilization Transition (VST).
  - {{contrast}}.StatWithVst.accepted.csv: all the gene differential analyses
    passed [the independent filtering](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilttheory) of DESeq2.

