"""Snakemake script for SnpSift with dbSNP"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake.utils import makedirs

# STDout contains the vcf file, it shall not be captured
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Extra parameters default value is an empty string
extra = snakemake.params.get("extra", "")

shell(
    "SnpSift annotate "                  # Tool and its subcommand
    "{extra} "                           # Extra parameters
    "{snakemake.input.database} "        # Path to phastCons database
    "{snakemake.input.calls} "           # Path to input vcf
    "> {snakemake.output.calls} "        # Path to annotated vcf
    "{log}"                              # Logging behaviour
)
