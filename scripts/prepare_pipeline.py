#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

"""
This script aims to prepare both design.tsv and config.yaml processed by
the vcf-annotate-SnpEff-SnpSift pipeline

For the design TSV file, it iterates over the given directory and lists
all VCF files. Their paths are saved in the dedicated TSV file.

For the configuration yaml file, it parses command line arguments and
saves them in the correct yaml format.

You can test this script with:
pytest -vv prepare_pipeline.py

Usage example:
python3.8 prepare_pipeline.py

IGR Human example:
python3.8 prepare_pipeline.py
"""


import argparse
import logging
import os
import os.path
import pandas
import pytest
import shlex
import sys
import yaml

from pathlib import Path
from typing import Any, Dict, Generator, List, Optional


# Building custom class for help formatter
class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """


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
        "VCF_DIR",
        help="Path to the directory containing vcf files",
        type=str
    )

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
        "-i", "--indexes",
        help="Also search for VCF indexes",
        action="store_true"
    )

    main_parser.add_argument(
        "-r", "--recursive",
        help="Recursively search in sub-directories for fastq files",
        action="store_true"
    )

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

    main_parser.add_argument(
        "--force",
        help="Over-write output file if they exist",
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
def parse_args(args: Any) -> argparse.ArgumentParser:
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


@pytest.mark.parametrize(
    "command, expected", [
        ("vcf GWASCat.tsv GenesSets.gmt dbNSFP.tsv --indexes",
         argparse.Namespace(
             GWASCat='GWASCat.tsv',
             GeneSets='GenesSets.gmt',
             cold_storage=[' '],
             dbNSFP='dbNSFP.tsv',
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
             workdir='.',
             recursive=False,
             indexes=True,
             VCF_DIR='vcf',
             force=False
        )),
        ("vcf GWASCat.tsv GenesSets.gmt dbNSFP.tsv",
         argparse.Namespace(
             GWASCat='GWASCat.tsv',
             GeneSets='GenesSets.gmt',
             cold_storage=[' '],
             dbNSFP='dbNSFP.tsv',
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
             workdir='.',
             recursive=False,
             indexes=False,
             VCF_DIR='vcf',
             force=False
        ))
    ]
)
def test_parse_args(command: str,
                    expected: argparse.Namespace) -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    print(shlex.split(command))
    assert parse_args(shlex.split(command)) == expected


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
    """
    result_dict = {
        "design": os.path.join(os.path.abspath(args.workdir), "design.tsv"),
        "config": os.path.join(os.path.abspath(args.workdir), "config.yaml"),
        "workdir": os.path.abspath(args.workdir),
        "force": args.force,
        "design_params": {
            "index": args.indexes,
            "recursive": args.recursive,
            "vcf_dir": os.path.abspath(args.VCF_DIR)
        },
        "threads": args.threads,
        "singularity_docker_image": args.singularity,
        "cold_storage": args.cold_storage,
        "ref": {
            "GWASCat": os.path.abspath(args.GWASCat),
            "dbNSFP": os.path.abspath(args.dbNSFP),
            "GeneSets": os.path.abspath(args.GeneSets)
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
    logging.debug(result_dict)
    return result_dict


@pytest.mark.parametrize(
    "command, expected", [
        ("vcf GWASCat.tsv GenesSets.gmt dbNSFP.vcf --indexes",
         {
             "design": os.path.abspath("./design.tsv"),
             "config": os.path.abspath("./config.yaml"),
             "workdir": os.path.abspath("."),
             "threads": 1,
             "force": False,
             "design_params": {
                 "index": True,
                 "recursive": False,
                 "vcf_dir": os.path.abspath("vcf")
             },
             "singularity_docker_image": 'docker://continuumio/'
                                         'miniconda3:4.4.10',
             "cold_storage": [' '],
             "ref": {
                 "GWASCat": os.path.abspath("GWASCat.tsv"),
                 "dbNSFP": os.path.abspath("dbNSFP.vcf"),
                 "GeneSets": os.path.abspath("GenesSets.gmt")
             },
             "organism": "GRCh38.86",
             "workflow": {
                 "multiqc": True
             },
             "params": {
                 "snpeff_extra": "-v",
                 "snpsift_dbNSFP_extra": "-v",
                 "snpsift_GWASCat_extra": "-v",
                 "snpsift_GeneSets_extra": "-v"
             }
         })
    ]
)
def test_args_to_dict(command: str, expected: Dict[str, str]) -> None:
    """
    This function simply tests the args_to_dict function with expected output

    Example:
    >>> pytest -v prepare_config.py -k test_args_to_dict
    """
    assert args_to_dict(parse_args(shlex.split(command))) == expected


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


# Looking for vcf files
def search_vcf(vcf_dir: Path,
               recursive: bool = False,
               exts: List[str] = (".vcf", ".vcf.gz")) \
               -> Generator[str, str, None]:
    """
    Iterate over a directory and search for vcf files (or any file ending
    with given extension list)

    Parameters:
        vcf_dir     Path        Path to the vcf directory in which to search
        recursive   bool        A boolean, weather to search recursively in
                                sub-directories (True) or not (False)
        exts        List[str]   A list of extensions used to identify vcf (or
                                any other file)

    Return:
                    Generator[str, str, None]       A Generator of paths
    """
    for path in vcf_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_vcf(path, recursive, exts)
            else:
                continue

        if path.name.endswith(exts):
            yield path


@pytest.mark.parametrize(
    "vcfd, rec, ext, expected", [
        (Path("tests/vcfs/empty"), False, None,
         [Path("tests/vcfs/empty/test.vcf")]),
        (Path("tests/vcfs/empty"), False, ("tbi"), []),
        (Path("tests/vcfs"), True, None,
         [Path("tests/vcfs/empty/test.vcf")])
    ]
)
def test_search_vcfs(vcfd: str,
                     rec: bool,
                     ext: List[str],
                     expected: List[str]) -> None:
    """
    This function tests the above test_search_vcfs with variable arguments
    """
    if ext is not None:
        assert list(search_vcf(vcfd, rec, ext)) == expected
    else:
        assert list(search_vcf(vcfd, rec)) == expected



# Turning the VCF list into a dictionnary
def classify_vcf(vcf_files: List[Path],
                 index_files: Optional[List[Path]] = None) \
                 -> Dict[str, Dict[str, str]]:
    """
    Return a dictionnary with identified vcf files (indexed or not)

    Parameters:
        vcf_files   List[Path]      A list of paths to iterate over
        index_files List[Path]      A list of paths to iterate over

    Return:
                    Dict[str, Path] A dictionnary: for each Sample ID, the ID
                                    is repeated alongside with the vcf and its
                                    corresponding index
    """
    vcf_dict = None
    if (index_files is not None) and (len(vcf_files) == len(index_files)):
        logging.debug("Matching number of indexes and vcf files")
        vcf_dict = {
            vcf.name: {
                "Sample_id": vcf.stem,
                "VCF_File": vcf.absolute(),
                "VCF_Index": index.absolute()
            }
            for vcf, index in zip(sorted(vcf_files), sorted(index_files))
        }
    else:
        logging.debug("No VCF index taken into account")
        vcf_dict = {
            vcf.name: {
                "Sample_id": vcf.stem,
                "VCF_File": vcf.absolute()
            }
            for vcf in vcf_files
        }
    return vcf_dict


@pytest.mark.parametrize(
    "vcf_list, index_list, expected", [
        ([Path("/path/to/file.vcf")], [Path("/path/to/file.vcf.tbi")],
         {"file.vcf": {"Sample_id": "file",
                       "VCF_File": Path("/path/to/file.vcf"),
                       "VCF_Index": Path("/path/to/file.vcf.tbi")}}),
        ([Path("/path/to/file.vcf.gz")], [Path("/path/to/file.vcf.gz.tbi")],
         {"file.vcf.gz": {"Sample_id": "file.vcf",
                          "VCF_File": Path("/path/to/file.vcf.gz"),
                          "VCF_Index": Path("/path/to/file.vcf.gz.tbi")}}),
        ([Path("/path/to/file.vcf.gz")], None,
         {"file.vcf.gz": {"Sample_id": "file.vcf",
                          "VCF_File": Path("/path/to/file.vcf.gz")}}),
        ([Path("path/to/file.vcf.gz")], None,
         {"file.vcf.gz": {"Sample_id": "file.vcf",
                          "VCF_File": Path("path/to/file.vcf.gz").absolute()}}),
    ]
)
def test_classify_vcf(vcf_list: List[Path],
                      index_list: Optional[List[Path]],
                      expected: Dict[str, str]) -> None:
    """
    This function tests the above classify_vcf function with multiple
    parameters.
    """
    if index_list is None:
        assert classify_vcf(vcf_list) == expected
    else:
        assert classify_vcf(vcf_list, index_list) == expected


def build_config(args: argparse.ArgumentParser) -> Dict[str, str]:
    """
    Build and save config file, then return config itself.
    """
    logging.info("Preparing configuration file")
    config = args_to_dict(args)
    output_path = Path(config["config"])

    if output_path.exists() and args.force is False:
        raise FileExistsError("Output file already exists: {config['config']}")

    with output_path.open('w') as config_file:
        logging.debug(f"Saving results to {config['config']}")
        config_file.write(dict_to_yaml(config))

    return config


def build_design(config: Dict[str, str]) -> None:
    """
    Build and save design file
    """
    vcf_files = sorted(
        list(search_vcf(
            Path(config["design_params"]["vcf_dir"]),
            config["design_params"]["recursive"]
        ))
    )

    index_files = sorted(
        list(search_vcf(
            Path(config["design_params"]["vcf_dir"]),
            config["design_params"]["recursive"],
            ("tbi")
        ))
    )

    logging.debug("Head of alphabeticaly sorted list of vcf files:")
    logging.debug([str(i) for i in vcf_files[0:5]])
    logging.debug([str(i) for i in index_files[0:5]])

    # Building a dictionnary of vcf (+ index?) and identifiers
    vcf_dict = classify_vcf(vcf_files, index_files)

    # Using Pandas to handle TSV output (yes pretty harsh I know)
    data = pandas.DataFrame(vcf_dict).T
    logging.debug("Head of vcf list table:")
    logging.debug(data.head())

    design_path = Path(config["design"])
    if design_path.exists() and config["force"] is False:
        raise FileExistsError(f"The file {design_path} already exists.")

    logging.debug(f"Saving results to {design_path}")
    data.to_csv(design_path, sep="\t", index=False)


# Running programm if script not imported
if __name__ == '__main__':
    # Pasing command line
    args = parse_args(sys.argv[1:])

    # Build logging object and behaviour
    os.makedirs("logs/prepare", exist_ok=True)
    logging.basicConfig(
        filename="logs/prepare/pipeline.log",
        filemode="w",
        level=logging.DEBUG
    )

    try:
        config = build_config(args)
        build_design(config)
    except Exception as e:
        logging.exception("%s", e)
        raise

    logging.info("Process over")
    sys.exit(0)
