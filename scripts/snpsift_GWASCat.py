"""Snakemake script for SnpSift GWAS Catalog"""

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

# In SnpSift GWASCat, one should give either -g parameter or -db parameter
# Not both.
db = snakemake.input.GWASCat is None
g_in_extra = " -g " not in extra

# n any of the previous cases ...
if db == g_in_extra:
    raise ValueError(
        "You should specify either a genome ('-g <genome>' in extra)"
        " **OR** give a database path (GWAS Catalog in input)"
    )

# Optional path to GWASCat database
try:
    db_arg = f"-db {snakemake.input.GWASCat}"
except KeyError:
    db_arg = ""

shell(
    "SnpSift gwasCat "                # Tool and its subcommand
    "{extra} "                        # Extra parameters
    "{db_arg} "                       # Path to GWAS Catalog
    "{snakemake.input.calls} "        # Path to input vcf
    "> {snakemake.output.calls} "     # Path to annotated vcf
    "{log}"                           # Logging behaviour
)
