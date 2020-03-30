#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This script contains functions that are to be called by any other scripts in
this pipeline.
"""

import argparse    # Argument parsing
import os          # Handle global OS operations
import os.path     # Handle paths and file system
import pandas      # Handle large datasets and tsv files
import pytest      # Testing behaviour

from typing import Any, Dict


def vcf_link(design: pandas.DataFrame) -> Dict[str, str]:
    """
    This function takes the "samples" described in config and returns
    a dictionnary with:
    sample file name : sample path
    """
    return {
        os.path.basename(vcf): os.path.realpath(vcf)
        for vcf in design["VCF_File"]
    }


@pytest.mark.parametrize(
    "design, expected", [
        (pandas.DataFrame({"VCF_File": {"S1": "/path/to/file.vcf"}}),
         {"file.vcf": "/path/to/file.vcf"}),
        (pandas.DataFrame({"VCF_File": {"S1": "path/to/file.vcf"}}),
         {"file.vcf": os.path.abspath("path/to/file.vcf")}),
        (pandas.DataFrame({"VCF_File": {"S1": "/path/to/file.vcf",
                                        "S2": "/other/file2.vcf"}}),
         {"file.vcf": "/path/to/file.vcf", "file2.vcf": "/other/file2.vcf"}),
        (pandas.DataFrame({"VCF_File": {"S1": "path/to/file.vcf.gz"}}),
         {"file.vcf.gz": os.path.abspath("path/to/file.vcf.gz")}),
    ]
)
def test_vcf_link(design: pandas.DataFrame, expected: Dict[str, str]) -> None:
    """
    Test the function vcf_link right above, with multiple parameters
    """
    print(design, vcf_link(design))
    assert vcf_link(design) == expected


def vcf_root(design: pandas.DataFrame) -> Dict[str, str]:
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
                base = os.path.basename(vcf)[:-(len(ext) + 1)]
                result[base] = f"raw_data/{os.path.basename(vcf)}"
                break
        else:
            raise ValueError(f"Could not remove ext: {vcf}")

    return result


@pytest.mark.parametrize(
    "design, expected", [
        (pandas.DataFrame({"VCF_File": {"S1": "/path/to/file.vcf"}}),
         {"file": "raw_data/file.vcf"}),
        (pandas.DataFrame({"VCF_File": {"S1": "path/to/file.vcf"}}),
         {"file": "raw_data/file.vcf"}),
        (pandas.DataFrame({"VCF_File": {"S1": "/path/to/file.vcf",
                                        "S2": "/other/file2.vcf"}}),
         {"file": "raw_data/file.vcf", "file2": "raw_data/file2.vcf"}),
        (pandas.DataFrame({"VCF_File": {"S1": "path/to/file.vcf.gz"}}),
         {"file": "raw_data/file.vcf.gz"}),
    ]
)
def test_vcf_root(design: pandas.DataFrame, expected: Dict[str, str]) -> None:
    """
    This function tests the above vcf_root funciton with multiple parameters
    """
    print(design, vcf_root(design))
    assert vcf_root(design) == expected


def named_refs_paths(config: Dict[str, Any]) -> Dict[str, str]:
    """
    This function returns the name of the databases and it used link
    within the pipeline
    """
    result = {}
    if (dbNSFP := config["ref"].get("dbNSFP", None)) is not None:
        result[os.path.basename(dbNSFP)] = os.sep.join(
            ["databases", os.path.basename(dbNSFP)]
        )
    else:
        raise FileNotFoundError("Could not find dbNSFP vcf file.")

    tbi = os.path.realpath(dbNSFP) + ".tbi"
    result[os.path.basename(tbi)] = os.sep.join(
        ["databases", os.path.basename(tbi)]
    )

    if (GWASCat := config["ref"].get("GWASCat", None)) is not None:
        result[os.path.basename(GWASCat)] = os.sep.join(
            ["databases", os.path.basename(GWASCat)]
        )
    else:
        raise FileNotFoundError("A GWAS Catalog tsv file is required.")

    if (GeneSets := config["ref"].get("GeneSets", None)) is not None:
        result[os.path.basename(GeneSets)] = os.sep.join(
            ["databases", os.path.basename(GeneSets)]
        )
    else:
        raise FileNotFoundError("Could not find GeneSets GMT file.")

    return result


@pytest.mark.parametrize(
    "conf, xpecd", [
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf": "databases/dbNSFP.vcf",
          "dbNSFP.vcf.tbi": "databases/dbNSFP.vcf.tbi",
          "GWASCat.tsv": "databases/GWASCat.tsv",
          "GeneSets.gmt": "databases/GeneSets.gmt"}),
        ({"ref": {"dbNSFP": "path/to/dbNSFP.vcf", "GWASCat": "path/to/GWASCat.tsv",
          "GeneSets": "path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf": "databases/dbNSFP.vcf",
          "dbNSFP.vcf.tbi": "databases/dbNSFP.vcf.tbi",
          "GWASCat.tsv": "databases/GWASCat.tsv",
          "GeneSets.gmt": "databases/GeneSets.gmt"}),
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf.gz", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf.gz": "databases/dbNSFP.vcf.gz",
          "dbNSFP.vcf.gz.tbi": "databases/dbNSFP.vcf.gz.tbi",
          "GWASCat.tsv": "databases/GWASCat.tsv",
          "GeneSets.gmt": "databases/GeneSets.gmt"})
    ]
)
def test_named_refs_paths(conf: Dict[str, Any], xpecd: Dict[str, str]) -> None:
    """
    This function tests the above named_refs_paths with multiple arguments
    """
    assert named_refs_paths(conf) == xpecd


def id_refs_paths(config: Dict[str, Any]) -> Dict[str, str]:
    """
    This function returns the id of the databases and it used link
    within the pipeline
    """
    result = {}
    if (dbNSFP := config["ref"].get("dbNSFP", None)) is not None:
        result["dbNSFP"] = os.sep.join(
            ["databases", os.path.basename(dbNSFP)]
        )
    else:
        raise FileNotFoundError("Could not find dbNSFP vcf file.")

    tbi = os.path.realpath(dbNSFP) + ".tbi"
    result["dbNSFP_tbi"] = os.sep.join(
        ["databases", os.path.basename(tbi)]
    )

    if (GWASCat := config["ref"].get("GWASCat", None)) is not None:
        result["GWASCat"] = os.sep.join(
            ["databases", os.path.basename(GWASCat)]
        )
    else:
        raise FileNotFoundError("A GWAS Catalog tsv file is required.")

    if (GeneSets := config["ref"].get("GeneSets", None)) is not None:
        result["GeneSets"] = os.sep.join(
            ["databases", os.path.basename(GeneSets)]
        )
    else:
        raise FileNotFoundError("Could not find GeneSets GMT file.")

    return result


@pytest.mark.parametrize(
    "conf, xpecd", [
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP": "databases/dbNSFP.vcf",
          "dbNSFP_tbi": "databases/dbNSFP.vcf.tbi",
          "GWASCat": "databases/GWASCat.tsv",
          "GeneSets": "databases/GeneSets.gmt"}),
        ({"ref": {"dbNSFP": "path/to/dbNSFP.vcf", "GWASCat": "path/to/GWASCat.tsv",
          "GeneSets": "path/to/GeneSets.gmt"}},
         {"dbNSFP": "databases/dbNSFP.vcf",
          "dbNSFP_tbi": "databases/dbNSFP.vcf.tbi",
          "GWASCat": "databases/GWASCat.tsv",
          "GeneSets": "databases/GeneSets.gmt"}),
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf.gz", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP": "databases/dbNSFP.vcf.gz",
          "dbNSFP_tbi": "databases/dbNSFP.vcf.gz.tbi",
          "GWASCat": "databases/GWASCat.tsv",
          "GeneSets": "databases/GeneSets.gmt"})
    ]
)
def test_id_refs_paths(conf: Dict[str, Any], xpecd: Dict[str, str]) -> None:
    """
    This function tests the above named_refs_paths with multiple arguments
    """
    assert id_refs_paths(conf) == xpecd


def link_refs_paths(config: Dict[str, Any]) -> Dict[str, str]:
    """
    This function returns the name of the databases and it realpath
    within the pipeline
    """
    result = {}
    if (dbNSFP := config["ref"].get("dbNSFP", None)) is not None:
        result[os.path.basename(dbNSFP)] = os.path.realpath(dbNSFP)
    else:
        raise FileNotFoundError("Could not find dbNSFP vcf file.")

    tbi = os.path.realpath(dbNSFP) + ".tbi"
    result[os.path.basename(tbi)] = tbi

    if (GWASCat := config["ref"].get("GWASCat", None)) is not None:
        result[os.path.basename(GWASCat)] = os.path.realpath(GWASCat)
    else:
        raise FileNotFoundError("A GWAS Catalog tsv file is required.")

    if (GeneSets := config["ref"].get("GeneSets", None)) is not None:
        result[os.path.basename(GeneSets)] = os.path.realpath(GeneSets)
    else:
        raise FileNotFoundError("Could not find GeneSets GMT file.")

    return result


@pytest.mark.parametrize(
    "conf, xpecd", [
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf": "/path/to/dbNSFP.vcf",
          "dbNSFP.vcf.tbi": "/path/to/dbNSFP.vcf.tbi",
          "GWASCat.tsv": "/path/to/GWASCat.tsv",
          "GeneSets.gmt": "/path/to/GeneSets.gmt"}),
        ({"ref": {"dbNSFP": "path/to/dbNSFP.vcf", "GWASCat": "path/to/GWASCat.tsv",
          "GeneSets": "path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf": os.path.abspath("path/to/dbNSFP.vcf"),
          "dbNSFP.vcf.tbi": os.path.abspath("path/to/dbNSFP.vcf.tbi"),
          "GWASCat.tsv": os.path.abspath("path/to/GWASCat.tsv"),
          "GeneSets.gmt": os.path.abspath("path/to/GeneSets.gmt")}),
        ({"ref": {"dbNSFP": "/path/to/dbNSFP.vcf.gz", "GWASCat": "/path/to/GWASCat.tsv",
          "GeneSets": "/path/to/GeneSets.gmt"}},
         {"dbNSFP.vcf.gz": "/path/to/dbNSFP.vcf.gz",
          "dbNSFP.vcf.gz.tbi": "/path/to/dbNSFP.vcf.gz.tbi",
          "GWASCat.tsv": "/path/to/GWASCat.tsv",
          "GeneSets.gmt": "/path/to/GeneSets.gmt"})
    ]
)
def test_link_refs_paths(conf: Dict[str, Any], xpecd: Dict[str, str]) -> None:
    """
    This function tests the above named_refs_paths with multiple arguments
    """
    assert link_refs_paths(conf) == xpecd
