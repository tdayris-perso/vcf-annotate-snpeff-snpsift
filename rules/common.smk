"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""


from typing import Dict, List, Tuple  # Type hinting for developpers
from snakemake.utils import validate  # Validate config and design formats

import os.path as op    # Path and file system manipulation
import os               # OS related operations
import pandas           # Deal with TSV files (design)
import sys              # System related operations


try:
    from common_vass import *
except ImportError:
    print(locals())
    raise

# Snakemake-Wrappers version
swv = "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/0.50.4"
# github prefix
git = "https://raw.githubusercontent.com/tdayris-perso/snakemake-wrappers"

# Loading configuration
if config == dict():
    configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pandas.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")

report: "../report/general.rst"


def get_vcf_w(wildcards) -> str:
    """
    Based on sample given as wildcars, this function returns vcf path
    """
    return vcf_root_dict[wildcards]


def get_targets(multiqc: bool = False) -> Dict[str, str]:
    """
    This function returns a dictionnary of all expected files at the end
    of this pipeline
    """
    targets = {}

    if (config["workflow"]["multiqc"] is True) and (multiqc is False):
        targets["multiqc"] = "qc/multiqc_report.html"

    if multiqc is False:
        targets["snpsift_GeneSets"] = expand(
            "snpsift/GeneSets/{sample}.vcf.gz",
            sample=sample_id_list
        )

        targets["snpeff_html"] = expand(
            "snpeff/report/{sample}.html",
            sample=sample_id_list
        )

        if config["workflow"]["multiqc"] is True:
            targets["multiqc"] = "qc/multiqc_report.html"

    targets["snpeff_csv"] = expand(
        "snpeff/stats/{sample}.csv",
        sample=sample_id_list
    )

    return targets


sample_id_list = design["Sample_id"].tolist()
vcf_link_dict = vcf_link(design)
vcf_root_dict = vcf_root(design)
ref_link_dict = link_refs_paths(config)
ref_pack_dict = named_refs_paths(config)
refs_dict = id_refs_paths(config)
targets_dict = get_targets()
