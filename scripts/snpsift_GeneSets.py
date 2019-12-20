"""Snakemake script for SnpSift Genes Sets from MSigDB"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.utils import makedirs

# STDout contains the vcf file, it shall not be captured
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extra parameters default value is an empty string
extra = snakemake.params.get("extra", "")

shell(
    "SnpSift geneSets "                  # Tool and its subcommand
    "{extra} "                           # Extra parameters
    "{snakemake.input.geneSets} "        # Path to MSigDB genes sets
    "{snakemake.input.calls} "           # Path to input vcf
    "> {snakemake.output.calls} "        # Path to annotated vcf
    "{log}"                              # Logging behaviour
)
