Material and Methods:
#####################

Annotation was performed by two tools from the same team: SnpEff, and SnpSift. No information was deleted from initial VCF files, no filters were applied. Quality reports were aggregated with MultiQC. The whole pipeline was powered by Snakemake, and the Snakemake-Wrappers project.

* SnpEff optional arguments: `{{snakemake.config.params.snpeff}}`
* SnpSift dnSNFP optional arguments: `{{snakemake.config.params.snpsift_dbSNFP}}`
* SnpSift GWASCat optional arguments: `{{snakemake.config.params.snpsift_GWASCat}}`
* SnpSift GeneSets optional arguments: `{{snakemake.config.params.snpsift_GeneSets}}`

If you need any other information, please read the `Frequently Asked questions <https://github.com/tdayris-perso/vcf-annotate-snpeff-snpsift#frequently-asked-questions-by-my-fellow-biologists-on-this-pipeline>`_ , then contact your bioinformatician if you're still in trouble.

Citations:
##########

SnpEff
  Cingolani, Pablo, et al. "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3." Fly 6.2 (2012): 80-92.

  http://snpeff.sourceforge.net/SnpEff.html

SnpSift
  Ruden, Douglas Mark, et al. "Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift." Frontiers in genetics 3 (2012): 35.

  http://snpeff.sourceforge.net/SnpSift.version_4_0.html

MultiQC
  EWELS, Philip, MAGNUSSON, Måns, LUNDIN, Sverker, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016, vol. 32, no 19, p. 3047-3048.

  https://multiqc.info/

Snakemake
  Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

  https://snakemake.readthedocs.io/
  https://snakemake-wrappers.readthedocs.io/
