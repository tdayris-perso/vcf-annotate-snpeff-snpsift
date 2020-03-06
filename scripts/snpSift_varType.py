"""Snakemake script for SnpSift varType subcommand"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.utils import makedirs

# STDout contains the vcf file, it shall not be captured
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "SnpSift varType "                  # Tool and its subcommand
    "{snakemake.input.calls} "           # Path to input vcf
    "> {snakemake.output.calls} "        # Path to annotated vcf
    "{log}"                              # Logging behaviour
)
