# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os

import numpy as np


from pathlib import Path
from collections import defaultdict

from do_alignment import align_reads, parse_alignments


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


def get_denominators(results):
    """Get denominators to normalise abundances.

    Args:
        results (dict): results

    Returns:
        numpy.ndarray: sum of columns in the metrics aka denominators for normalisation

    """
    return np.sum([results[element] for element in results], axis=0)


def add_relative_abundance(level_results, normalizer):
    # add relative abundance next to raw count for each of the file(s) in the analysis
    for level in level_results:
        relative_abundance = list((level_results[level]/normalizer) * 100)
        level_results[level] = list(level_results[level]) + relative_abundance

    return level_results


def aggregate_level(results, position, normalizer):
    """Aggregate abundance of subsystem level and add relative abundance.

    Args:
        results (dict): Path to results
        position (int): Position of level in the results
        normalizer (numpy.ndarray): normalizer denominators

    Returns:
        dict: Aggregated result targeting chosen subsystem level

    """
    level_results = defaultdict(list)

    for all_levels in results:
        level = all_levels.split("\t")[position]
        abundance = results[all_levels]
        level_results[level].append(abundance)

    level_results = {temp_level: np.sum(level_results[temp_level], axis=0) for temp_level in level_results}

    return add_relative_abundance(level_results, normalizer)


def get_subsystems(translation_file):
    """Create lookup table from primary key to subsystems levels 1, 2, and 3.

    Args:
        translation_file (PosixPath): Path to file with subsystems information

    Returns:
        dict: lookup table primary key to subsystem levels

    """
    subsystems_translation = {}
    with open(translation_file) as database_file:
        database_reader = csv.reader(database_file, delimiter='\t')
        next(database_reader, None)
        for row in database_reader:
            subsystems_translation[row[0]] = "\t".join(row[1:])

    return subsystems_translation


def write_results(results, header, output_name, query_path, database):
    """Write results in tabular format.

    Args:
        results (dict): dict with results to be written
        header (list): header to be wrritten
        output_name (str): Path to output
        query_path (str): Path to query
        database (str): Database used

    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        # run info
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow([""])

        # subsystem and files header
        writer.writerow(header)
        for row in sorted(results):
            if sum(results[row]) > 0:
                writer.writerow(row.split("\t") + results[row])


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
    minimum_identity = float(args.minimum_identity)
    minimum_alignment = int(args.minimum_alignment)
    threads = args.threads
    evalue = args.evalue
    database = args.database.split("_")[-1]
    amino_acid = args.amino_acid
    fast_mode = args.fast

    # other metrics
    keep_alignment = args.keep_alignment
    normalise_output = int(args.normalise_output)
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
        results = defaultdict(list)

        subsystems_translation = get_subsystems(Path(WORK_DIRECTORY,"db/database_PKs.txt"))

        # get fasta/fastq files
        query_files = is_wanted_file([temp_query for temp_query in os.listdir(queries_folder)])
        for counter, temp_query in enumerate(query_files):
            LOGGER.info("1.{}) Working on: {}".format(counter + 1, temp_query))
            LOGGER.info("   Aligning sequences in {} to {} using {}".format(temp_query, database, aligner))
            alignment_name = align_reads(Path(queries_folder, temp_query), output_directory, aligner, database, WORK_DIRECTORY)
            LOGGER.info("   Parsing Alignments")
            sample_position = query_files.index(temp_query)
            results = parse_alignments(alignment_name, results, normalise_output, len(query_files), sample_position,
                                        minimum_identity, minimum_alignment, subsystems_translation)

        # write results
        normalizer = get_denominators(results)
        header_files = query_files + ["{} %".format(x) for x in query_files]
        LOGGER.info('Writting results at {}'.format(output_directory))

        # write results for each of the levels
        for level in [1, 2, 3]:
            LOGGER.info('  Working on subsystem level {}'.format(level))
            temp_header = ["Subsystem {}".format(level)] + header_files

            temp_results = aggregate_level(results, level - 1, normalizer)
            output_file = "{}/{}.xls".format(output_directory, level)

            write_results(temp_results, temp_header, output_file, queries_folder, database)

        # write result for all the levels in one file
        LOGGER.info ('  Working on Combined output')
        temp_header = ["Subsystem Level 1", "Subsystem Level 2", "Subsystem Level 3"] + header_files
        output_file = "{}/combined.xls".format (output_directory)
        temp_results = add_relative_abundance(results, normalizer)
        write_results(temp_results, temp_header, output_file, queries_folder, database)


    LOGGER.info('Done'.format(output_directory))


if __name__ == "__main__":
    main()
