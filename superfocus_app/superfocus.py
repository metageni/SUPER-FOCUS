# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os
import random

from pathlib import Path
from collections import defaultdict

from do_alignment import align_reads


LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'
logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)
LOGGER = logging.getLogger(__name__)

WORK_DIRECTORY = 'superfocus_app'



def which(program_name):
    """python implementation of unix 'which' function.

    Args:
        program_name (str): Program name

    Returns:
        str: Program path

    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program_name)
    if fpath:
        if is_exe(program_name):
            return program_name
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program_name)
            if is_exe(exe_file):
                return exe_file

    return None


def is_wanted_file(queries):
    """Remove files from query files that not have extension .fasta/.fastq/.fna

    Args:
        queries (list): List with query names

    Returns:
        list: sorted list with only .fasta/.fastq/.fna files

    """
    queries = [query for query in queries if query.split(".")[-1].lower() in ["fna", "fasta", "fastq"]]
    queries.sort()

    return queries



def main():
    parser = argparse.ArgumentParser(description="SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data",
                                     epilog="python superfocus.py -q input_folder -dir output_dir")
    # basic parameters
    parser.add_argument("-q", "--query", help="Path to FAST(A/Q) file or directory with these files", required=True)
    parser.add_argument("-dir", "--output_directory",  help="Path to output files", required=True)
    parser.add_argument("-o", "--output_prefix",  help="Output prefix (Default: output)", default="output")

    # aligner related
    parser.add_argument("-a", "--aligner",  help="aligner choice (rapsearch, blast, diamond; default rapsearch)", default="rapsearch")
    parser.add_argument("-mi", "--minimum_identity",  help="minimum identity (default 60 perc)", default="60")
    parser.add_argument("-ml", "--minimum_alignment",  help="minimum alignment (amino acids) (default: 15)", default="15")
    parser.add_argument("-t", "--threads",  help="Number Threads used in the k-mer counting (Default: 4)", default="4")
    parser.add_argument("-e", "--evalue",  help="e-value (default 0.00001)", default="0.00001")
    parser.add_argument("-db", "--database",  help="database (DB_90, DB_95, DB_98, or DB_100; default DB_98)", default="DB_90")
    parser.add_argument("-p", "--amino_acid",  help="amino acid input; 0 nucleotides; 1 amino acids (default 0)", default="0.00001")
    parser.add_argument("-f", "--fast",  help="runs RAPSearch2 or DIAMOND on fast mode - 0 (False) / 1 (True) (default: 1)", default="1")

    # extra
    parser.add_argument("-k", "--keep_alignment",  help="keep original tabular output. 0 delete it / 1 keep it (default 0)", default="0")
    parser.add_argument("-n", "--normalise_output",  help="normalises each query counts based on number of hits; 0 doesn't normalize; 1 normalizes (default: 1)", default="1")
    parser.add_argument("-m", "--focus",  help="runs FOCUS; 1 does run; 0 does not run: default 0", default="0")
    parser.add_argument("-r", "--annotation",  help="use only the subsystems in the organisms predicted by -focus ncbi / rast annotation  (default: ncbi)", default="ncbi")
    parser.add_argument("-d", "--work_directory",  help="Work directory (Default: superfocus_app)", default="superfocus_app")
    parser.add_argument("-b", "--alternate_directory",  help="Alternate directory for your databases", default="")

    args = parser.parse_args()

    # basic parameters
    queries_folder = Path(args.query)
    prefix = args.output_prefix
    output_directory = Path(args.output_directory)

    # alignment related
    aligner = args.aligner.lower()
    minimum_identity = args.minimum_identity
    minimum_alignment = args.minimum_alignment
    threads = args.threads
    evalue = args.evalue
    database = args.database.split("_")[-1]
    amino_acid = args.amino_acid
    fast_mode = args.fast

    # other metrics
    keep_alignment = args.keep_alignment
    normalise_output = args.normalise_output
    run_focus = args.focus
    annotation = args.annotation
    WORK_DIRECTORY = Path(args.alternate_directory) if args.alternate_directory else Path(args.work_directory)

    LOGGER.info("SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data")

    # check if output_directory is exists - if not, creates it
    if not output_directory.exists():
        Path(output_directory).mkdir(parents=True, mode=511)
        LOGGER.info("OUTPUT: {} does not exist - just created it :)".format(output_directory))

    # check if at least one of the queries is valid
    if not queries_folder.is_dir():
        LOGGER.critical("QUERY: {} is not a directory".format(queries_folder))

    # check if at least one of the queries is valid
    elif is_wanted_file(os.listdir(queries_folder)) == []:
        LOGGER.critical("QUERY: {} does not have any Fasta/Fna/Fastq file".format(queries_folder))

    # check if aligner is valid
    elif aligner not in ["diamond", "rapsearch", "blast"]:
        LOGGER.critical("ALIGNER: {} is not a valid aligner. Please chose among (diamond, rapsearch, or blast)".format(aligner))

    # check if aligner exists
    elif not which(aligner):
        LOGGER.critical("ALIGNER: {} is not in the path of your system".format(aligner))

    # check if query is exists
    elif not queries_folder.exists():
        LOGGER.critical("QUERY: {} does not exist".format(queries_folder))

    # check if work directory exists
    elif WORK_DIRECTORY != WORK_DIRECTORY or not WORK_DIRECTORY.exists():
        LOGGER.critical("WORK_DIRECTORY: {} does not exist".format(WORK_DIRECTORY))

    else:
        # get fasta/fastq files
        query_files = is_wanted_file([temp_query for temp_query in os.listdir(queries_folder)])
        for counter, temp_query in enumerate(query_files):
            LOGGER.info("1.{}) Working on: {}".format(counter + 1, temp_query))
            LOGGER.info("   Aligning sequences in {} to {} using {}".format(temp_query, database, aligner))
            align_reads(Path(queries_folder, temp_query), output_directory, aligner, database, WORK_DIRECTORY)
            LOGGER.info("   Parsing Alignments")

    LOGGER.info('Done'.format(output_directory))


if __name__ == "__main__":
    main()
