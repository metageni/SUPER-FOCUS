#!/usr/bin/env python3
"""Pure stateless utility functions for SUPER-FOCUS."""

import csv
import logging
import os
import sys
import itertools

import numpy as np

from pathlib import Path
from collections import defaultdict


logger = logging.getLogger(__name__)


def is_valid_number(value):
    """Check if input is a valid >= 0 int or float.

    Args:
        value (str): Value to be checked.

    Returns:
        bool: True if valid >= 0 number else False.
    """
    try:
        return float(value) >= 0
    except (ValueError, TypeError):
        return False


def is_wanted_file(queries):
    """Filter query list to only .fasta/.fastq/.fna files, sorted.

    Args:
        queries (list): List of Path objects.

    Returns:
        list: Sorted list with only .fasta/.fastq/.fna files.
    """
    queries = [q for q in queries if q.name.split(".")[-1].lower() in ("fna", "fasta", "fastq")]
    queries.sort()
    return queries


def get_denominators(results):
    """Get column-wise sums of results for relative abundance normalisation.

    Args:
        results (dict): Results dict mapping function keys to count lists.

    Returns:
        numpy.ndarray: Sum of columns (denominators for normalisation).
    """
    return np.sum([results[element] for element in results], axis=0)


def add_relative_abundance(level_results, normalizer):
    """Append relative abundance columns to each result row.

    Args:
        level_results (dict): Results to be updated.
        normalizer (numpy.ndarray): Per-sample total counts.

    Returns:
        dict: Results with relative abundance appended.
    """
    for level in level_results:
        relative_abundance = np.divide(list(level_results[level]), normalizer, where=normalizer != 0)
        relative_abundance *= 100
        level_results[level] = list(level_results[level]) + list(relative_abundance)
    return level_results


def aggregate_level(results, position, normalizer):
    """Aggregate counts at a given subsystem level and append relative abundance.

    Args:
        results (dict): Full results dict.
        position (int): Tab-delimited field index to aggregate on.
        normalizer (numpy.ndarray): Per-sample total counts.

    Returns:
        dict: Aggregated results with relative abundance.
    """
    level_results = defaultdict(list)
    for all_levels in results:
        level = all_levels.split("\t")[position]
        level_results[level].append(results[all_levels])
    level_results = {k: np.sum(v, axis=0) for k, v in level_results.items()}
    return add_relative_abundance(level_results, normalizer)


def get_subsystems(translation_file):
    """Load subsystem lookup table from database_PKs.txt.

    Args:
        translation_file (Path): Path to tab-delimited PKs file.

    Returns:
        dict: Mapping of primary key to tab-joined subsystem levels.
    """
    subsystems_translation = {}
    with open(translation_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)
        for row in reader:
            subsystems_translation[row[0]] = "\t".join(row[1:])
    return subsystems_translation


def write_results(results, header, output_name, query_path, database, aligner):
    """Write functional profile results to a tab-delimited file.

    Args:
        results (dict): Results dict mapping function keys to count lists.
        header (list): Column header row.
        output_name (str): Output file path.
        query_path (str): Query path string for metadata header.
        database (str): Database name for metadata header.
        aligner (str): Aligner name for metadata header.

    Returns:
        None
    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow(["Aligner used: {}".format(aligner)])
        writer.writerow([""])
        writer.writerow(header)
        for row in sorted(results):
            if sum(results[row]) > 0:
                writer.writerow(row.split("\t") + list(map(str, results[row])))


def write_binning(binning_result, output_name, query_path, database, aligner):
    """Write per-read binning results to a tab-delimited file.

    Args:
        binning_result (defaultdict): Nested dict of query→read→hit list.
        output_name (str): Output file path.
        query_path (str): Query path string for metadata header.
        database (str): Database name for metadata header.
        aligner (str): Aligner name for metadata header.

    Returns:
        None
    """
    with open(output_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        writer.writerow(["Query: {}".format(query_path)])
        writer.writerow(["Database used: {}".format(database)])
        writer.writerow(["Aligner used: {}".format(aligner)])
        writer.writerow([""])
        writer.writerow(["Sample name", "Read Name",
                         "Subsystem Level 1", "Subsystem Level 2", "Subsystem Level 3", "Function",
                         "Identity %", "Alignment Length", "E-value"])
        for query_name in binning_result:
            for read_name in binning_result[query_name]:
                temp_row = binning_result[query_name][read_name]
                temp_row = list(temp_row for temp_row, _ in itertools.groupby(temp_row))
                for row_temp in temp_row:
                    writer.writerow([query_name, read_name] + row_temp[-1].split("\t") + row_temp[:-1])


def subsample_reads(input_file, output_directory, number_of_reads):
    """Subsample the first N reads from a fasta/fastq file.

    Args:
        input_file (Path): Path to input fasta/q file.
        output_directory (str): Directory for the subsampled output file.
        number_of_reads (int): Number of reads to keep.

    Returns:
        Path: Path to the subsampled output file.
    """
    tmpoutput = Path(os.path.join(output_directory, input_file.name))
    if tmpoutput.exists():
        logger.critical(
            "%s already exists in %s — refusing to overwrite during subsampling",
            input_file.name, output_directory,
        )
        sys.exit(-1)

    logger.info("   Subsampling %s, selecting %d → %s", input_file, number_of_reads, tmpoutput)
    with open(input_file, 'r') as fin, open(tmpoutput, 'w') as out:
        seqcounter = 0
        if 'fna' in input_file.suffixes or '.fasta' in input_file.suffixes:
            for r in fin:
                if r.startswith('>'):
                    seqcounter += 1
                if seqcounter <= number_of_reads:
                    out.write(r)
        else:
            linecounter = 0
            for r in fin:
                if linecounter == 0:
                    seqcounter += 1
                linecounter += 1
                if linecounter == 4:
                    linecounter = 0
                if seqcounter <= number_of_reads:
                    out.write(r)
    return tmpoutput
