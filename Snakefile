import snakemake.utils  # Load snakemake API
import sys              # System related operations

# Python 3.7 is required
if sys.version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later.")

# Snakemake 5.4.2 at least is required
snakemake.utils.min_version("5.13.0")

include: "rules/common.smk"
include: "rules/copy.smk"
include: "rules/multiqc.smk"
include: "rules/snpeff.smk"
include: "rules/snpsift.smk"
include: "rules/vcf.smk"

workdir: config["workdir"]
containers: config["singularity_docker_image"]
localrules: copy_fastq, copy_extra

rule all:
    input:
        multiqc = "qc/multiqc_report.html",
        snpsift_GeneSets = expand(
            "snpsift/GeneSets/{sample}.vcf.gz",
            sample=sample_id_list
        ),
        snpsift_GeneSets_tbi = expand(
            "snpsift/GeneSets/{sample}.vcf.gz.tbi",
            sample=sample_id_list
        )
    message:
        "Finishing the VCF annotation pipeline"
