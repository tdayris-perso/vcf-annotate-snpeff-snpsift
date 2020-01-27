Material and Methods:
#####################

Annotation was performed by two tools from the same team: SnpEff, and SnpSift. No information was deleted from initial VCF files, no filters were applied. Quality reports were aggregated with MultiQC. The whole pipeline was powered by `Snakemake <https://snakemake.readthedocs.io/
https://snakemake-wrappers.readthedocs.io/>`_ , and the `Snakemake-Wrappers <https://snakemake.readthedocs.io/
https://snakemake-wrappers.readthedocs.io/>`_ project.

* SnpEff optional arguments: `{{snakemake.config.params.snpeff}}`
* SnpSift dnSNFP optional arguments: `{{snakemake.config.params.snpsift_dbSNFP}}`
* SnpSift GWASCat optional arguments: `{{snakemake.config.params.snpsift_GWASCat}}`
* SnpSift GeneSets optional arguments: `{{snakemake.config.params.snpsift_GeneSets}}`

If you need any other information, please read the `Frequently Asked questions <https://github.com/tdayris-perso/vcf-annotate-snpeff-snpsift#frequently-asked-questions-by-my-fellow-biologists-on-this-pipeline>`_ , then contact your bioinformatician if you're still in trouble.

Citations:
##########

This pipeline stands on best practices found in multiple high impact papers, published in Nature, Cell, Bioinformatics, and others.

SnpEff
  Cingolani, Pablo, et al. "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3." Fly 6.2 (2012): 80-92.

  Why SnpEff? SnpEff is a very well known, powerfull and complete tool for annotating vairants. It has been cited more than 3500 times. It performs faster than Annovar, and remains completely open.

  http://snpeff.sourceforge.net/SnpEff.html

SnpSift
  Ruden, Douglas Mark, et al. "Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift." Frontiers in genetics 3 (2012): 35.

  SnpSift makes the pair with SnpEff. From the very same team, these tools almost never come without the other. SnpSift takes up-to-date annotations, on a wider range of annotations than Annovar, and remains completely open.

  http://snpeff.sourceforge.net/SnpSift.version_4_0.html

MultiQC
  EWELS, Philip, MAGNUSSON, Måns, LUNDIN, Sverker, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016, vol. 32, no 19, p. 3047-3048.

  Why MultiQC? MultiQC is a very efficient tool when it comes to quality gathering. It has been cited more than 500 times in a very wide range of journals including Nature, Bioinformatics, Cell, etc.

  https://multiqc.info/

Snakemake
  Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

  Why Snakemake? Snakemake is a very popular workflow manager in data science and bioinformatics. It has about three new citations per week within the scopes of biology, medicine and bioinformatics.

  https://snakemake.readthedocs.io/
  https://snakemake-wrappers.readthedocs.io/
