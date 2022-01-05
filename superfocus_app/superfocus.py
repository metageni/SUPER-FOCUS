#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import itertools
import logging
import os
import shutil
import sys
import tempfile

import numpy as np

from pathlib import Path
from shutil import which
from collections import defaultdict

from superfocus_app.do_alignment import align_reads, parse_alignments
from superfocus_app import version

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'


def is_wanted_file(queries):
    """Remove files from query files that not have extension .fasta/.fastq/.fna

    Args:
        queries (list): List with query names.

    Returns:
        list: Sorted list with only .fasta/.fastq/.fna files.

    """
    queries = [query for query in queries if query.name.split(".")[-1].lower() in ["fna", "fasta", "fastq"]]
    queries.sort()

    return queries


def get_denominators(results):
    """Get denominators to normalise abundances.

    Args:
        results (dict): Results.

    Returns:
        numpy.ndarray: Sum of columns in the metrics aka denominators for normalisation.

    """
    return np.sum([results[element] for element in results], axis=0)


def add_relative_abundance(level_results, normalizer):
    """Add relative abundance to results.

    Args:
        level_results (dict): Results to be updated.
        normalizer (numpy.ndarray): Normalizer denominators.

    Returns:
        dict: Results with relative abundance.

    """
    # add relative abundance next to raw count for each of the file(s) in the analysis
    for level in level_results:
        relative_abundance = np.divide(list(level_results[level]), normalizer, where=normalizer != 0)
        relative_abundance *= 100
        level_results[level] = list(level_results[level]) + list(relative_abundance)

    return level_results


def aggregate_level(results, position, normalizer):
    """Aggregate abundance of subsystem level and add relative abundance.

    Args:
        results (dict): Results.
        position (int): Position of level in the results.
        normalizer (numpy.ndarray): Normalizer denominators.

    Returns:
        dict: Aggregated result targeting chosen subsystem level.

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
        translation_file (PosixPath): Path to file with subsystems information.

    Returns:
        dict: Lookup table primary key to subsystem levels.

    """
    subsystems_translation = {}
    with open(translation_file) as database_file:
        database_reader = csv.reader(database_file, delimiter='\t')
        next(database_reader, None)
        for row in database_reader:
            subsystems_translation[row[0]] = "\t".join(row[1:])

    return subsystems_translation


def write_results(results, header, output_name, query_path, database, aligner):
    """Write results in tabular format.

    Args:
        results (collections.defaultdict): Dict with results to be written.
        header (list): Header to be written.
        output_name (str): Path to output.
        query_path (str): Path to query.
        database (str): Database used.
        aligner (str): Aligner name.

    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        # run info
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow(["Aligner used: {}".format(aligner)])
        writer.writerow([""])

        # subsystem and files header
        writer.writerow(header)
        for row in sorted(results):
            if sum(results[row]) > 0:
                writer.writerow(row.split("\t") + list(map(str, results[row])))


def write_binning(binning_result, output_name, query_path, database, aligner):
    """Write binning results in tabular format.

    Args:
        binning_result (collections.defaultdict): Dict with results to be written.
        output_name (str): Path to output.
        query_path (str): Path to query.
        database (str): Database used.
        aligner (str): Aligner name.

    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        # run info
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow(["Aligner used: {}".format(aligner)])
        writer.writerow([""])

        writer.writerow(["Sample name", "Read Name",
                         "Subsystem Level 1", "Subsystem Level 2", "Subsystem Level 3", "Function",
                         "Identity %", "Alignment Length", "E-value"])
        for query_name in binning_result:
            for read_name in binning_result[query_name]:
                # remove duplicates from list
                temp_row = binning_result[query_name][read_name]
                temp_row = list(temp_row for temp_row, _ in itertools.groupby(temp_row))
                for row_temp in temp_row:
                    row = [query_name, read_name] + row_temp[-1].split("\t") + row_temp[:-1]
                    writer.writerow(row)


def is_valid_number(value):
    """ Check if input if a valid >= 0 int or float.

    Args:
        value (str): Value to be checked.

    Returns:
        bool: True if valid >= 0 number else False.

    """
    try:
        if float(value) >= 0:
            return True
        else:
            return False
    except:
        return False


def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: Parsed arguments.

    """
    parser = argparse.ArgumentParser(description="SUPER-FOCUS: A tool for agile functional analysis of shotgun "
                                                 "metagenomic data.",
                                     epilog="superfocus -q input_folder -dir output_dir")
    parser.add_argument('-v', '--version', action='version', version='SUPER-FOCUS {}'.format(version))
    # basic parameters
    parser.add_argument("-q", "--query", help="Path to FAST(A/Q) file or directory with these files.", required=True,
                        action='append')
    parser.add_argument("-dir", "--output_directory", help="Path to output files", required=True)
    parser.add_argument("-o", "--output_prefix", help="Output prefix (Default: output).", default="output_")
    parser.add_argument("-tmp", "--temp_directory", help="specify an alternate temporary directory to use")

    # aligner related
    parser.add_argument("-a", "--aligner", help="aligner choice (rapsearch, diamond, blast, or mmseqs2; default rapsearch).",
                        default="rapsearch")
    parser.add_argument("-mi", "--minimum_identity", help="minimum identity (default 60 perc).", default="60")
    parser.add_argument("-ml", "--minimum_alignment", help="minimum alignment (amino acids) (default: 15).",
                        default="15")
    parser.add_argument("-t", "--threads", help="Number Threads used in the k-mer counting (Default: 4).",
                        default="4")
    parser.add_argument("-e", "--evalue", help="e-value (default 0.00001).", default="0.00001")
    parser.add_argument("-db", "--database", help="database (DB_90, DB_95, DB_98, or DB_100; default DB_90)",
                        default="DB_90")
    parser.add_argument("-p", "--amino_acid", help="amino acid input; 0 nucleotides; 1 amino acids (default 0).",
                        default="0")
    parser.add_argument("-f", "--fast", help="runs RAPSearch2 or DIAMOND on fast mode - 0 (False) / 1 (True) "
                        "(default: 1).", default="1")

    # extra
    parser.add_argument("-n", "--normalise_output", help="normalises each query counts based on number of hits; "
                                                         "0 doesn't normalize; 1 normalizes (default: 1).", default="1")
    parser.add_argument("-m", "--focus", help="runs FOCUS; 1 does run; 0 does not run: default 0.", default="0")
    parser.add_argument("-b", "--alternate_directory", help="Alternate directory for your databases.", default="")
    parser.add_argument('-d', '--delete_alignments', help='Delete alignments', action='store_true', required=False)
    parser.add_argument('-w', '--latency_wait', help='Add a delay (in seconds) between writing the file and reading it',
                        type=int, default=0)
    parser.add_argument('-l', '--log', help='Path to log file (Default: STDOUT).', required=False)

    return parser.parse_args()


def main():
    args = parse_args()

    # basic parameters
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
    del_alignments = args.delete_alignments

    # rename mmseqs2 to mmseqs. This is the name of the command, but the program is mmseqs!
    if aligner == 'mmseqs2':
        aligner = 'mmseqs'

    # other metrics
    normalise_output = int(args.normalise_output)
    run_focus = args.focus
    if args.alternate_directory:
        WORK_DIRECTORY = Path(args.alternate_directory)
    elif 'SUPERFOCUS_DB' in os.environ:
        WORK_DIRECTORY = Path(os.environ['SUPERFOCUS_DB'])
    else:
        WORK_DIRECTORY =  Path(__file__).parents[0]


    if args.log:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO, filename=args.log)
    else:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)

    logger = logging.getLogger(__name__)


    logger.info("SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data")

    # check if output_directory is exists - if not, creates it
    if not output_directory.exists():
        Path(output_directory).mkdir(parents=True, mode=511)
        logger.info("OUTPUT: {} does not exist - just created it :)".format(output_directory))

    # parse directory and/or query files
    query_files = []
    for f in args.query:
        p = Path(f)
        if p.is_dir():
            query_files += [Path(p, x) for x in os.listdir(p)]
        elif p.is_file():
            query_files.append(p)
    query_files = is_wanted_file(query_files)
    if query_files == []:
        logger.critical("QUERY: {} does not have any fasta/fna/fastq files".format(args.query))
        sys.exit(1)

    # find a temp directory location
    tmp = "/tmp"
    if args.temp_directory:
        tmp = args.temp_directory
        if not os.path.exists(tmp):
            logger.info(f"Creating temporary path: {tmp}")
            os.makedirs(tmp, exist_ok=True)
    elif 'TMPDIR' in os.environ:
        tmp = os.environ['TMPDIR']
        if not os.path.exists(tmp):
            logger.info(f"Creating temporary path: {tmp}")
            os.makedirs(tmp, exist_ok=True)
    else:
        sys.stderr.write(f"WARNING: Using {tmp} as the base temporary directory")
    tmpdir = tempfile.mkdtemp(dir=tmp)
    os.makedirs(tmpdir, exist_ok=True)
    logger.info(f"Using {tmpdir} as the temporary directory")

    # check if we can run focus
    if run_focus != '0':
        logger.critical("FOCUS: Running FOCUS is not avaliable on this version. "
                        "Please see https://github.com/metageni/FOCUS on how to run it")

    # check if amino_acid is valid
    elif aligner == 'blast' and amino_acid not in ['0', '1']:
        logger.critical("AMINO ACID OPTION: {} is not valid for --amino_acid. Only 0 or 1".format(amino_acid))

    # check if at database choice is valid
    elif database not in ["90", "95", "98", "100"]:
        logger.critical("DATABASE: DB_{} not valid. Choose DB_90/95/98/or 100".format(database))

    # check if aligner is valid
    elif aligner not in ["diamond", "rapsearch", "blast", "mmseqs"]:
        logger.critical("ALIGNER: {} is not a valid aligner. Please choose among (diamond, blast, rapsearch, or mmseqs2)".
                        format(aligner))

    # check if aligner exists
    elif not which(aligner) and aligner.lower() != "blast":
        logger.critical("ALIGNER: {} is not in the path of your system".format(aligner))

    # check if work directory exists
    elif WORK_DIRECTORY != WORK_DIRECTORY or not WORK_DIRECTORY.exists():
        logger.critical("WORK_DIRECTORY: {} does not exist".format(WORK_DIRECTORY))

    # check if number of threads are valid
    elif threads != "all" and not is_valid_number(threads):
        logger.critical("THREADS: {} is not a valid number of threads".format(threads))

    # check if evalue is valid
    elif not is_valid_number(evalue):
        logger.critical("E-VALUE: {} is not a valid evalue".format(evalue))

    else:
        results = defaultdict(list)
        binning_reads = defaultdict(lambda: defaultdict(list))

        subsystems_translation = get_subsystems(Path(WORK_DIRECTORY, "db/database_PKs.txt"))

        for counter, temp_query in enumerate(query_files):
            logger.info("1.{}) Working on: {}".format(counter + 1, temp_query))
            logger.info("   Aligning sequences in {} to {} using {}".format(temp_query, database, aligner))
            alignment_name = align_reads(temp_query,
                                         output_dir=output_directory, aligner=aligner,
                                         database=database, evalue=evalue,
                                         threads=threads, fast_mode=fast_mode,
                                         WORK_DIRECTORY=WORK_DIRECTORY,
                                         amino_acid=amino_acid, temp_folder=tmpdir,
                                         latency_delay=args.latency_wait)
            logger.info("   Parsing Alignments")
            sample_position = query_files.index(temp_query)
            results, binning_reads = parse_alignments(alignment_name, results, normalise_output, len(query_files),
                                                      sample_position, minimum_identity, minimum_alignment,
                                                      subsystems_translation, aligner, binning_reads, temp_query,
                                                      del_alignments)

        # write results
        normalizer = get_denominators(results)
        header_files = query_files + ["{} %".format(x) for x in query_files]
        logger.info('Writting results at {}'.format(output_directory))

        # write binning
        output_file = "{}/{}binning.xls".format(output_directory, prefix)
        write_binning(binning_reads, output_file, args.query, database, aligner)
        logger.info('  Working on writing binning')

        # write results for each of the levels
        for level in [1, 2, 3]:
            logger.info('  Working on subsystem level {}'.format(level))
            temp_header = ["Subsystem {}".format(level)] + header_files

            temp_results = aggregate_level(results, level - 1, normalizer)
            output_file = "{}/{}subsystem_level_{}.xls".format(output_directory, prefix, level)

            write_results(temp_results, temp_header, output_file, args.query, database, aligner)

        # write result for all the levels in one file
        logger.info('  Working on Combined output')
        temp_header = ["Subsystem Level 1", "Subsystem Level 2", "Subsystem Level 3", "Function"] + header_files
        output_file = "{}/{}all_levels_and_function.xls".format(output_directory, prefix)
        temp_results = add_relative_abundance(results, normalizer)
        write_results(temp_results, temp_header, output_file, args.query, database, aligner)

    # clean up our mess
    shutil.rmtree(tmpdir)

    logger.info('Done')


if __name__ == "__main__":
    main()
