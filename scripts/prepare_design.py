#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script aims to prepare the list of files to be processed
by the vcf-annotate-snpeff-snpsift pipeline

It iterates over a given directory, lists all vcf files,
and saves their paths in a TSV file.

You can test this script with:
pytest -v ./prepare_design.py

Usage example:
python3.7 ./prepare_design.py tests/reads

# Search in sub-directories:
python3.7 ./prepare_design.py tests --recursive
"""

import argparse           # Parse command line
import logging            # Traces and loggings
import logging.handlers   # Logging behaviour
import os                 # OS related activities
import pandas as pd       # Parse TSV files
import pytest             # Unit testing
import shlex              # Lexical analysis
import sys                # System related methods

from pathlib import Path                        # Paths related methods
from typing import Any, Dict, Generator, List, Optional  # Type hints

from common_script_vass import *

logger = setup_logging(logger="prepare_design.py")


# Processing functions
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

    Example:
    # Basic usage
    >>> search_vcf(Path("tests/vcfs/"))
    <generator object search_fq at 0xXXXXXXXXXXXX>

    # Recover list of files
    >>> list(search_vcf(Path("tests/", True)))
    [Path('tests/vcf/test.vcf')]

    # Recusively search for vcf indexes
    >>> search_vcf(Path("tests/vcfs"), recursive = True, exts = ["tbi"])
    """
    for path in vcf_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_vcf(path, recursive, exts)
            else:
                continue

        if path.name.endswith(exts):
            yield path


# Testing search_vcf
def test_search_vcfs() -> None:
    """
    This function tests the ability of the function "search_vcf" to find the
    vcf files in the given directory

    Example:
    pytest -v prepare_design.py -k test_search_vcf
    """
    path = Path("tests/vcfs/")
    expected = [
        Path("tests/vcfs/test.vcf")
    ]
    assert list(search_vcf(path)) == sorted(expected)


def test_search_vcfs_index() -> None:
    """
    This function tests the ability of the function "search_vcf" to find the
    vcf index files in the given directory

    Example:
    pytest -v prepare_design.py -k test_search_vcf
    """
    path = Path("tests/vcfs/")
    expected = []

    assert list(search_vcf(path, exts=("tbi"))) == sorted(expected)


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

    Example:
    # VCF with index example
    >>> classify_vcf([Path("file1.vcf"), Path("file2.vcf")],
                     [Path("file1.vcf.tbi"), Path("File2.vcf.tbi")])
    {'file1': {'VCF_File': Path("/path/to/file1.vcf"),
     'Sample_id': 'file1',
     'VCF_Index': Path('/path/to/file1.vcf.tbi')}

    # VCF without index example
    >>> classify_vcf([Path("file1.vcf"), Path("file2.vcf")])
    {'file1': {'VCF_File': Path("/path/to/file1.vcf"),
     'Sample_id': 'file1'}

    # Example with non-matching number of VCF and indexes
    >>> classify_vcf([Path("file1.vcf"), Path("file2.vcf")],
                     [Path("file1.vcf.tbi"), Path("File2.vcf.tbi")])
    {'file1': {'VCF_File': Path("/path/to/file1.vcf"),
     'Sample_id': 'file1'}
    """
    vcf_dict = None
    if (index_files is not None) and (len(vcf_files) == len(index_files)):
        logger.debug("Matching number of indexes and vcf files")
        vcf_dict = {
            vcf.name: {
                "Sample_id": vcf.stem,
                "VCF_File": vcf.absolute(),
                "VCF_Index": index.absolute()
            }
            for vcf, index in zip(sorted(vcf_files), sorted(index_files))
        }
    else:
        logger.debug("No VCF index taken into account")
        vcf_dict = {
            vcf.name: {
                "Sample_id": vcf.stem,
                "VCF_File": vcf.absolute()
            }
            for vcf in vcf_files
        }
    return vcf_dict


def test_classify_vcf_with_index() -> None:
    """
    This function takes input from the pytest decorator
    to test the classify_vcf function

    Example:
    pytest -v ./prepare_design.py -k test_classify_vcf_with_index
    """
    prefix = Path(__file__).parent.parent
    vcf_list = [prefix / "file1.vcf"]
    index_list = [prefix / "file1.vcf.tbi"]

    expected = {
        'file1.vcf': {
            'VCF_File': Path(prefix / "file1.vcf"),
            'Sample_id': 'file1',
            'VCF_Index': Path(prefix / 'file1.vcf.tbi')
        }
    }
    assert classify_vcf(vcf_list, index_list) == expected


def test_classify_vcf_no_index() -> None:
    """
    This function takes input from the pytest decorator
    to test the classify_vcf function

    Example:
    pytest -v ./prepare_design.py -k test_classify_vcf_no_index
    """
    prefix = Path(__file__).parent.parent
    vcf_list = [prefix / "file1.vcf"]
    index_list = []

    expected = {
        'file1.vcf': {
            'VCF_File': Path(prefix / "file1.vcf"),
            'Sample_id': 'file1'
        }
    }
    assert classify_vcf(vcf_list, index_list) == expected


def test_classify_vcf_non_matching_index() -> None:
    """
    This function takes input from the pytest decorator
    to test the classify_vcf function

    Example:
    pytest -v ./prepare_design.py -k test_classify_vcf_non_matching_index
    """
    prefix = Path(__file__).parent.parent
    vcf_list = [prefix / "file1.vcf", prefix / "file2.vcf"]
    index_list = [prefix / "file2.vcf.tbi"]

    expected = {
        'file1.vcf': {
            'VCF_File': Path(prefix / "file1.vcf"),
            'Sample_id': 'file1'
        },
        'file2.vcf': {
            "VCF_File": Path(prefix / "file2.vcf"),
            "Sample_id": "file2"
        }
    }
    assert classify_vcf(vcf_list, index_list) == expected


# Parsing command line arguments
# This function won't be tested
def parse_args(args: Any = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    Build a command line parser object

    Parameters:
        args    Any                 Command line arguments

    Return:
                ArgumentParser      Parsed command line object

    Example:
    >>> parse_args(shlex.split("/path/to/vcf --indexes"))
    Namespace(debug=False, output='design.tsv', path='/path/to/vcf',
    quiet=False, recursive=False, indexes=True)
    """
    # Defining command line options
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does not perform any magic. Check the result."
    )

    # Required arguments
    main_parser.add_argument(
        "path",
        help="Path to the directory containing vcf files",
        type=str
    )

    # Optional arguments
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
        "-o", "--output",
        help="Path to output file (default: %(default)s)",
        type=str,
        default="design.tsv"
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

    # Parsing command lines
    return main_parser.parse_args(args)


def test_parse_args() -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    options = parse_args(shlex.split("/path/to/vcf/dir/ --index"))
    expected = argparse.Namespace(
        debug=False,
        output='design.tsv',
        path='/path/to/vcf/dir/',
        quiet=False,
        recursive=False,
        indexes=True
    )
    assert options == expected


# Main function, the core of this script
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole preparation sequence

    Parameters:
        args    ArgumentParser      The parsed command line

    Example:
    >>> main(parse_args(shlex.split("/path/to/fasta/dir/")))
    """
    # Searching for vcf files and sorting them alphabetically
    vcf_files = sorted(
        list(search_vcf(Path(args.path), args.recursive))
    )
    index_files = sorted(
        list(search_vcf(Path(args.path), args.recursive, ("tbi")))
    )

    logger.debug("Head of alphabeticaly sorted list of vcf files:")
    logger.debug([str(i) for i in vcf_files[0:5]])
    logger.debug([str(i) for i in index_files[0:5]])

    # Building a dictionnary of vcf (+ index?) and identifiers
    vcf_dict = classify_vcf(vcf_files, index_files)

    # Using Pandas to handle TSV output (yes pretty harsh I know)
    data = pd.DataFrame(vcf_dict).T
    logger.debug("\n{}".format(data.head()))
    logger.debug("Saving results to {}".format(args.output))
    data.to_csv(args.output, sep="\t", index=False)


# Running programm if not imported
if __name__ == '__main__':
    # Parsing command line
    args = parse_args()

    try:
        logger.debug("Preparing design")
        main(args)
    except Exception as e:
        logger.exception("%s", e)
        sys.exit(1)
    sys.exit(0)
