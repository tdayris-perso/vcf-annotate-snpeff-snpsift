#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script aims to prepare the configuration file used
by the vcf-annotate-snpeff-snpsift pipeline

It goes through the arguments passed in command line and
builds a yaml formatted text file used as a configuration
file for the snakemake pipeline.

You can test this script with:
pytest -v ./prepare_config.py

Usage example:
# Whole pipeline
python3.7 ./prepare_config.py /path/to/GWASCat.tsv \
/path/to/GenesSets.gmt \
/path/to/dbNSFP.tsv

# Whole pipeline, verbose mode activated
python3.7 ./prepare_config.py /path/to/GWASCat.tsv \
/path/to/GenesSets.gmt \
/path/to/dbNSFP.tsv \
--verbose
"""


import argparse             # Parse command line
import logging              # Traces and loggings
import logging.handlers     # Logging behaviour
import os                   # OS related activities
import pytest               # Unit testing
import shlex                # Lexical analysis
import sys                  # System related methods
import yaml                 # Parse Yaml files

from pathlib import Path             # Paths related methods
from typing import Dict, Any         # Typing hints

from common import *

logger = setup_logging(logger="prepare_config.py")


def parser() -> argparse.ArgumentParser:
    """
    Build the argument parser object
    """
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does not make any magic. Please check the prepared"
               " configuration file!"
    )

    # Parsing positional argument
    main_parser.add_argument(
        "GWASCat",
        help="Path to the GWAS Catalog file",
        type=str,
    )

    main_parser.add_argument(
        "GeneSets",
        help="Path to the GeneSets GMT-formatted file",
        type=str,
    )

    main_parser.add_argument(
        "dbNSFP",
        help="Path to the dbNSFP TSV-formatted file",
        type=str,
    )

    # Parsing optional arguments
    main_parser.add_argument(
        "--design",
        help="Path to design file (default: %(default)s)",
        type=str,
        metavar="PATH",
        default="design.tsv"
    )

    main_parser.add_argument(
        "--workdir",
        help="Path to working directory (default: %(default)s)",
        type=str,
        metavar="PATH",
        default="."
    )

    main_parser.add_argument(
        "--threads",
        help="Maximum number of threads used (default: %(default)s)",
        type=int,
        default=1
    )

    main_parser.add_argument(
        "--singularity",
        help="Docker/Singularity image (default: %(default)s)",
        type=str,
        default="docker://continuumio/miniconda3:4.4.10"
    )

    main_parser.add_argument(
        "--cold-storage",
        help="Space separated list of absolute path to "
             "cold storage mount points (default: %(default)s)",
        nargs="+",
        type=str,
        default=[" "]
    )

    main_parser.add_argument(
        "--snpeff-extra",
        help="Extra parameters for SnpEff step (default: %(default)s)",
        type=str,
        default="-v"
    )

    main_parser.add_argument(
        "--organism",
        help="Organism used by SnpEff",
        type=str,
        default="GRCh38.86"
    )

    main_parser.add_argument(
        "--snpsift-GWASCat-extra",
        help="Extra parameters for SnpSift GWASCat (default: %(default)s)",
        type=str,
        default="-v"
    )

    main_parser.add_argument(
        "--snpsift-dbNSFP-extra",
        help="Extra parameters for SnpSift dbNSFP (default: %(default)s)",
        type=str,
        default="-v"
    )

    main_parser.add_argument(
        "--snpsift-GeneSets-extra",
        help="Extra parameters for SnpSift GeneSets (default: %(default)s)",
        type=str,
        default="-v"
    )

    main_parser.add_argument(
        "--no-multiqc",
        help="Disable multiqc aggregation in this pipeline",
        action="store_true",
        default=False
    )

    # Logging options
    log = main_parser.add_mutually_exclusive_group()
    log.add_argument(
        "-d", "--debug",
        help="Set logging in debug mode",
        default=False,
        action='store_true'
    )

    log.add_argument(
        "-q", "--quiet",
        help="Turn off logging behaviour",
        default=False,
        action='store_true'
    )

    return main_parser


# Argument parsing functions
def parse_args(args: Any = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    This function parses command line arguments

    Parameters
        args     Any             All command line arguments

    Return
                ArgumentParser   A object designed to parse the command line

    Example:
    >>> parse_args(shlex.split("/path/to/GWASCat.tsv /path/to/GenesSets.gmt
    /path/to/dnNSFP.tsv"))
    """
    return parser().parse_args(args)


def test_parse_args() -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    options = parse_args(shlex.split("/path/to/GWASCat.tsv "
                                     "/path/to/GenesSets.gmt "
                                     "/path/to/dbNSFP.tsv"))
    expected = argparse.Namespace(
        GWASCat='/path/to/GWASCat.tsv',
        GeneSets='/path/to/GenesSets.gmt',
        cold_storage=[' '],
        dbNSFP='/path/to/dbNSFP.tsv',
        debug=False, design='design.tsv',
        organism="GRCh38.86",
        no_multiqc=False,
        quiet=False,
        singularity='docker://continuumio/miniconda3:4.4.10',
        snpeff_extra='-v',
        snpsift_GWASCat_extra='-v',
        snpsift_GeneSets_extra='-v',
        snpsift_dbNSFP_extra='-v',
        threads=1,
        workdir='.'
    )
    assert options == expected


# Building pipeline configuration from command line
def args_to_dict(args: argparse.ArgumentParser) -> Dict[str, Any]:
    """
    Parse command line arguments and return a dictionnary ready to be
    dumped into yaml

    Parameters:
        args        ArgumentParser      Parsed arguments from command line

    Return:
                    Dict[str, Any]      A dictionnary containing the parameters
                                        for the pipeline

    Examples:
    >>> example_options = parse_args("/path/to/fasta")
    >>> args_to_dict(example_options)
    {
        'design': '/path/to/design',
        'workdir': '/path/to/workdir',
        'threads': 100,
        'singularity_docker_image': 'singularity_image',
        'cold_storage': ['/path/cold/one', '/path/cold/two'],
        'ref': {
            'GWASCat': '/path/to/GWASCat.tsv',
            'dbNSFP': '/path/to/dbNSFP.tsv',
            'GeneSets': '/path/to/GenesSets.gmt'
        },
        'organism': "GRCh38.86",
        'workflow': {'multiqc': True},
        'params': {
            'snpeff_extra': '-v',
            'snpsift_dbNSFP_extra': '-v',
            'snpsift_GWASCat_extra': '-v',
            'snpsift_GeneSets_extra': '-v'
        }
    }
    """
    result_dict = {
        "design": args.design,
        "workdir": args.workdir,
        "threads": args.threads,
        "singularity_docker_image": args.singularity,
        "cold_storage": args.cold_storage,
        "ref": {
            "GWASCat": args.GWASCat,
            "dbNSFP": args.dbNSFP,
            "GeneSets": args.GeneSets
        },
        "organism": args.organism,
        "workflow": {
            "multiqc": not args.no_multiqc
        },
        "params": {
            "snpeff_extra": args.snpeff_extra,
            "snpsift_dbNSFP_extra": args.snpsift_dbNSFP_extra,
            "snpsift_GWASCat_extra": args.snpsift_GWASCat_extra,
            "snpsift_GeneSets_extra": args.snpsift_GeneSets_extra
        }
    }
    logger.debug(result_dict)
    return result_dict


def test_args_to_dict() -> None:
    """
    This function simply tests the args_to_dict function with expected output

    Example:
    >>> pytest -v prepare_config.py -k test_args_to_dict
    """
    options = parse_args(shlex.split(
        "/path/to/GWASCat.tsv "
        "/path/to/GenesSets.gmt "
        "/path/to/dbNSFP.tsv "
        "--design /path/to/design "
        "--workdir /path/to/workdir "
        "--threads 100 "
        "--singularity singularity_image "
        "--cold-storage /path/cold/one /path/cold/two "
        "--debug "
    ))

    expected = {
        'design': '/path/to/design',
        'workdir': '/path/to/workdir',
        'threads': 100,
        'singularity_docker_image': 'singularity_image',
        'cold_storage': ['/path/cold/one', '/path/cold/two'],
        'ref': {
            'GWASCat': '/path/to/GWASCat.tsv',
            'dbNSFP': '/path/to/dbNSFP.tsv',
            'GeneSets': '/path/to/GenesSets.gmt'
        },
        'organism': "GRCh38.86",
        'workflow': {'multiqc': True},
        'params': {
            'snpeff_extra': '-v',
            'snpsift_dbNSFP_extra': '-v',
            'snpsift_GWASCat_extra': '-v',
            'snpsift_GeneSets_extra': '-v'
        }
    }
    assert args_to_dict(options) == expected


# Yaml formatting
def dict_to_yaml(indict: Dict[str, Any]) -> str:
    """
    This function makes the dictionnary to yaml formatted text

    Parameters:
        indict  Dict[str, Any]  The dictionnary containing the pipeline
                                parameters, extracted from command line

    Return:
                str             The yaml formatted string, directly built
                                from the input dictionnary

    Examples:
    >>> import yaml
    >>> example_dict = {
        "bar": "bar-value",
        "foo": ["foo-list-1", "foo-list-2"]
    }
    >>> dict_to_yaml(example_dict)
    'bar: bar-value\nfoo:\n- foo-list-1\n- foo-list-2\n'
    >>> print(dict_to_yaml(example_dict))
    bar: bar-value
    foo:
    - foo-list-1
    - foo-list-2
    """
    return yaml.dump(indict, default_flow_style=False)


def test_dict_to_yaml() -> None:
    """
    This function tests the dict_to_yaml function with pytest

    Example:
    >>> pytest -v prepare_config.py -k test_dict_to_yaml
    """
    expected = 'bar: bar-value\nfoo:\n- foo-list-1\n- foo-list-2\n'
    example_dict = {
        "bar": "bar-value",
        "foo": ["foo-list-1", "foo-list-2"]
    }
    assert dict_to_yaml(example_dict) == expected


# Core of this script
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole configuration sequence

    Parameters:
        args    ArgumentParser      The parsed command line

    Example:
    >>> main(parse_args(shlex.split("/path/to/fasta")))
    """
    # Building pipeline arguments
    logger.debug("Building configuration file:")
    config_params = args_to_dict(args)
    output_path = Path(args.workdir) / "config.yaml"

    # Saving as yaml
    with output_path.open("w") as config_yaml:
        logger.debug(f"Saving results to {str(output_path)}")
        config_yaml.write(dict_to_yaml(config_params))


# Running programm if not imported
if __name__ == '__main__':
    # Parsing command line
    args = parse_args()

    try:
        logger.debug("Preparing configuration")
        main(args)
    except Exception as e:
        logger.exception("%s", e)
        sys.exit(1)
    sys.exit(0)
