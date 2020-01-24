"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""

from snakemake.utils import validate
from typing import Any, Dict, List, Tuple

import os.path as op    # Path and file system manipulation
import os               # OS related operations
import pandas as pd     # Deal with TSV files (design)
import sys              # System related operations

# Snakemake-Wrappers version
swv = "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/0.49.0"
# github prefix
git = "https://raw.githubusercontent.com/tdayris-perso/snakemake-wrappers"

# Loading configuration
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pd.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")

report: "../report/general.rst"


def vcf_link() -> Dict[str, str]:
    """
    This function takes the "samples" described in config and returns
    a dictionnary with:
    sample file name : sample path
    """
    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    return {
        op.basename(vcf): op.realpath(vcf)
        for vcf in design["VCF_File"]
    }


def vcf_root() -> Dict[str, str]:
    """
    This function takes the VCF file list and returns the root
    name corresponding to a VCF file in the following dictionnary structure:
    sample name: sample link path
    """
    # For now, bz2 compression is not taken into account.
    possible_ext = ("vcf", "vcf.gz")

    # Build final result
    result = {}
    for vcf in design["VCF_File"]:
        # I always love writing these crazy for-break-else!
        for ext in possible_ext:
            if vcf.endswith(ext):
                # Extension removal
                base = op.basename(vcf)[:-(len(ext) + 1)]
                result[base] = f"raw_data/{op.basename(vcf)}"
                break
        else:
            raise ValueError(f"Could not remove ext: {vcf}")

    return result


def ref_pack() -> Tuple[Dict[str, str]]:
    """
    Returns three dictionaries with references paths
    """
    real = {}
    pack = {}
    named = {}

    dbNSFP = config["ref"].get("dbNSFP", None)
    if dbNSFP is not None:
        real[op.basename(dbNSFP)] = op.realpath(dbNSFP)
        pack[op.basename(dbNSFP)] = os.sep.join(
            ["databases", op.basename(dbNSFP)]
        )
        named["dbNSFP"] = os.sep.join(
            ["databases", op.basename(dbNSFP)]
        )
        tbi = op.realpath(dbNSFP) + ".tbi"
        if op.exists(tbi):
            real[op.basename(tbi)] = tbi
            pack[op.basename(tbi)] = os.sep.join(
                ["databases", op.basename(tbi)]
            )
            named["dbNSFP_tbi"] = os.sep.join(
                ["databases", op.basename(tbi)]
            )
        else:
            raise FileNotFoundError("dbNSFP index not found.")
    else:
        raise FileNotFoundError("A dbNSFP TSV file is required.")

    GWASCat = config["ref"].get("GWASCat", None)
    if GWASCat is not None:
        real[op.basename(GWASCat)] = op.realpath(GWASCat)
        pack[op.basename(GWASCat)] = os.sep.join(
            ["databases", op.basename(GWASCat)]
        )
        named["GWASCat"] = os.sep.join(
            ["databases", op.basename(GWASCat)]
        )
    else:
        raise FileNotFoundError("A GWAS catalog TSV file is required.")

    GeneSets = config["ref"].get("GeneSets", None)
    if GeneSets is not None:
        real[op.basename(GeneSets)] = op.realpath(GeneSets)
        pack[op.basename(GeneSets)] = os.sep.join(
            ["databases", op.basename(GeneSets)]
        )
        named["GeneSets"] = os.sep.join(
            ["databases", op.basename(GeneSets)]
        )
    else:
        raise FileNotFoundError("A GeneSets GMT file is required.")

    return real, pack, named


def sample_id() -> List[str]:
    """
    Return the list of samples identifiers
    """
    return design["Sample_id"].tolist()


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


sample_id_list = sample_id()
vcf_link_dict = vcf_link()
vcf_root_dict = vcf_root()
ref_link_dict, ref_pack_dict, refs_dict = ref_pack()
targets_dict = get_targets()
