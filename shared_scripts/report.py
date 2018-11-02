#!/usr/bin/env python3
from snakemake.utils import report
#print (snakemake.input)
report("""
A RNA-Seq workflow to map reads to a genome+transcriptome and to detect differitial expressed genes.
============================================================= 
Workflow:
-------------------------------------------------------------
1. Raw reads for each sample were concatenated into fastq file (only R1 single-end reads).
#. Reads were mapped to genome+transcriptome using STAR
#. subreads::featureCounts were used to count reads to annotated genes.
#. DESeq2 were used for differential expression.

A detailed presentation of the workflow:

.. image:: dag.png
    :target: dag.png

Result Files: 
-------------------------------------------------------------
- `config.yaml <config.yaml>`_: the configurantion file for the workflow.
  Including genome and transcriptome sources, contrasts for each DESeq2
  analysis, and parameters for other steps of the workflow.

- `sampleInfo.csv <sampleInfo.csv>`_: the sampleInfo to map sampleName to group

- hg38/

  - \*.bam: the alignment file for human genome.

- hg38.gene_count/

  - geneCount.csv: a gene count matrix of gene_name|gene_id.version x sampleName
  - geneRpkm.csv: a Reads Per Kbp per Million total mapped (RPKM) matrix.
  - `geneSummary.csv <hg38.gene_count/geneSummary.csv>`_: a summary of count by gene_type. 
  - `geneSummary.pdf <hg38.gene_count/geneSummary.pdf>`_: a summary of count by gene_type. 

- hg38.deseq2/: :raw-html:`<br>`
  The DESeq2 analyses of all the contrasts, non-paired. (see config.yaml)

  - plotPca_top500.pdf: the PCA plot using the top 500 genes with large variance.
  - expr.vst_blind.csv: the count matrix normalized with DESeq2
    Variance Stabilization Transition (VST).
  - {{contrast}}.StatWithVst.accepted.csv: all the gene differential analyses
    passed `the independent filtering <http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilttheory>`_ of DESeq2.

    - avgLog2Rpkm: the average of log2 RPKM for each sample in the contrast.
    - baseMean: the average reads of each sample after normalization with total reads.
    - log2FoldChange: log2(contrast2/contrast1)
    - stat: the statistics based on the distribution or count of reads
      (negative binormial)
    - pvalue: the pvalue base on the statistics above
    - padj: p value adjusted for multiple comparison using FDR method.
    - The last columns: the expr on log2 scale for each sample (normalized
      by DESeq2::VarianceStabilizationTransformation) 

  - {{contrast}}.StatWithVst.padj_05.addRpkm.csv: the gene differential
    analyses that is significant by padj<.05, with the VST expr
    values at the last columns, and the average log2 RPKM at the
    first column.

  - {{contrast}}.StatWithVst.padj_05.addRpkm.png: a plot to help decide
    filtering criteria.


Tips for downstream analyses:
-------------------------------------------------------------
- Cutoff criteria of the significant genes
  
  The default cutoff in the report is padj < .05.  You can add log2FoldChange
  and/or log2AvgRpkm to your criteria.

- Pathway/Geneset analyses:

  `DAVID <https://david.ncifcrf.gov/list.jsp>`_ can be a good start.  Just
  paste the genesymbols of your significant genes (up- and down- regulated as
  seperate lists).

- Comparison of significant genes from different contrasts

  Venndiagrams can be a good start, try `intervene on-line app
  <https://asntech.shinyapps.io/intervene/>`_ or `venny
  <http://bioinfogp.cnb.csic.es/tools/venny/>`_. Note that it is suggested to
  compare the upregulated and downregulated seperately.


.. - hg38.deseq2.paired_ic/: The DESeq2 analyses of ic samples, paired. (see
..   config.yaml)
.. - mm10/
.. - mm10.gene_count/
.. - mm10.deseq2.grouped_all/
.. - mm10.deseq2.paired_ic/
.. - mm10.deseq2.paired_wt/

""", snakemake.output[0]) #, **input)

