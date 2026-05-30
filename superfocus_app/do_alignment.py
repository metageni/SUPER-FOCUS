#!/usr/bin/env python3
"""Alignment execution and parsing functions for SUPER-FOCUS."""

import csv
import logging
import os
import sys
import subprocess

from pathlib import Path
from collections import defaultdict

_PID = os.getpid()

csv.field_size_limit(100000000)

logger = logging.getLogger(__name__)


def normalise_counts(data):
    """Normalise hit counts so they sum to 1.

    Args:
        data (dict): Mapping of function key to raw count.

    Returns:
        dict: Normalised counts summing to 1.0.
    """
    sum_hits = sum(data.values())
    return {level: data[level] / sum_hits for level in data}


def update_results(results, sample_index, data, normalise, number_samples):
    """Accumulate per-sample hit counts into the results dict.

    Args:
        results (defaultdict): Accumulated results across all samples.
        sample_index (int): Index of the current sample.
        data (dict): Hit counts for the current read group.
        normalise (int): 1 to normalise counts, 0 to keep raw.
        number_samples (int): Total number of samples.

    Returns:
        defaultdict: Updated results.
    """
    if data:
        if normalise and len(data) > 1:
            data = normalise_counts(data)
        for level in data:
            if level not in results:
                results[level] = [0] * number_samples
            results[level][sample_index] += data[level]
    return results


def align_reads(query, output_dir, aligner, database, evalue, threads, fast_mode,
                WORK_DIRECTORY, amino_acid, temp_folder, latency_delay=0):
    """Run an external aligner against the SUPER-FOCUS database.

    Args:
        query (Path): Path to query FAST(A/Q) file.
        output_dir (Path): Directory for alignment output.
        aligner (str): One of diamond, mmseqs, rapsearch, blast.
        database (str): Database cluster size (e.g. '90').
        evalue (str): E-value threshold.
        threads (str): Thread count or 'all'.
        fast_mode (str): '1' for fast mode, '0' for sensitive.
        WORK_DIRECTORY (Path): Root directory containing db/static/.
        amino_acid (str): '1' for amino acid input, '0' for nucleotide.
        temp_folder (str): Temporary directory for aligner scratch files.
        latency_delay (int): Seconds to wait after alignment (NFS latency).

    Returns:
        str: Path to the produced alignment file (.m8).
    """
    output_name = "{}/{}_{}_alignments".format(output_dir, query.parts[-1], _PID)
    database = "{}_clusters".format(database)
    blast_mode = 'blastp' if amino_acid == '1' else 'blastx'

    if aligner == "diamond":
        database_diamond = "{}/db/static/diamond/{}.db".format(WORK_DIRECTORY, database)
        cmd = [
            "diamond", blast_mode,
            "-d", database_diamond,
            "-q", query,
            "-o", f"{output_name}.m8",
            "-f", "6",
            "-t", temp_folder,
            "-p", threads,
            "-e", evalue,
        ]
        if fast_mode != "1":
            cmd.append("--sensitive")
        try:
            retcode = subprocess.call(cmd)
            if retcode != 0:
                logger.error("diamond exited with code %d", retcode)
                sys.exit(retcode)
        except OSError as e:
            logger.error("diamond execution failed: %s", e)
            sys.exit(1)
        output_name = "{}.m8".format(output_name)

    elif aligner == 'mmseqs':
        database_mmseqs = f"{WORK_DIRECTORY}/db/static/mmseqs2/{database}.db"
        output_name = f"{output_name}.m8"
        cmd = [
            "mmseqs", "easy-search",
            query, database_mmseqs, output_name, temp_folder,
            "--threads", threads,
            "-e", evalue,
        ]
        if fast_mode == "1":
            cmd += ["-s", "1.0"]
        try:
            retcode = subprocess.call(cmd)
            if retcode != 0:
                logger.error("mmseqs2 exited with code %d", retcode)
                sys.exit(retcode)
        except OSError as e:
            logger.error("mmseqs2 execution failed: %s", e)
            sys.exit(1)

    elif aligner == "rapsearch":
        mode_rapsearch = "T" if fast_mode == "1" else "F"
        database_rapsearch = "{}/db/static/rapsearch2/{}.db".format(WORK_DIRECTORY, database)
        os.system('rapsearch -a {} -q {} -d {} -o {} -v 250 -z {} -e {} -b 0 -s f'.format(
            mode_rapsearch, query, database_rapsearch, output_name, threads, evalue))
        output_name = "{}.m8".format(output_name)

    elif aligner == "blast":
        database_blast = "{}/db/static/blast/{}.db".format(WORK_DIRECTORY, database)
        os.system('{} -db {} -query {} -out {} -outfmt 6 -evalue {} -max_target_seqs 250 -num_threads {}'.format(
            blast_mode, database_blast, query, output_name, evalue, threads))

    else:
        logger.error("Unknown aligner: %s", aligner)
        sys.exit(1)

    return output_name


def parse_alignments(alignment, results, normalise, number_samples, sample_index,
                     minimum_identity, minimum_alignment, subsystems_translation,
                     aligner, binning_reads, query_name, delete_alignments):
    """Parse a tabular alignment file and accumulate functional hits.

    Args:
        alignment (str): Path to the alignment file (.m8 tabular format).
        results (defaultdict): Accumulated results across all samples.
        normalise (int): 1 to normalise per-read counts, 0 for raw.
        number_samples (int): Total number of samples.
        sample_index (int): Index of the current sample.
        minimum_identity (float): Minimum percent identity threshold.
        minimum_alignment (int): Minimum alignment length threshold.
        subsystems_translation (dict): PK → tab-joined subsystem levels.
        aligner (str): Aligner name (affects identity parsing for mmseqs).
        binning_reads (defaultdict): Per-read binning accumulator.
        query_name (str): Name of the query file (used as binning key).
        delete_alignments (bool): Delete the alignment file after parsing.

    Returns:
        tuple: (results, binning_reads) or defaultdict(int) if file is empty/missing.
    """
    if not os.path.exists(alignment) or os.stat(alignment).st_size == 0:
        return defaultdict(int)

    previous_read_name = None
    best_evalue = None

    with open(alignment, encoding='ISO-8859-1') as f:
        reader = csv.reader(f, delimiter='\t')
        if aligner == "rapsearch":
            [next(reader, None) for _ in range(5)]

        temp_results = defaultdict(int)
        for row in reader:
            current_read_name = row[0]
            temp_mi = float(row[2])
            if aligner == 'mmseqs' and 0 <= temp_mi <= 1:
                temp_mi *= 100
            temp_ml = float(row[3])
            current_evalue = row[10]

            current_hit = row[1].split("__")
            current_subsystem_id = current_hit[1]
            current_function_name = current_hit[-1].replace("\n", "").replace("\r", "")
            aggregate_levels = subsystems_translation[current_subsystem_id] + "\t" + current_function_name

            if current_read_name != previous_read_name:
                if previous_read_name:
                    update_results(results, sample_index, temp_results, normalise, number_samples)
                temp_results = defaultdict(int)
                best_evalue = current_evalue

            if temp_mi >= minimum_identity and temp_ml >= minimum_alignment and current_evalue == best_evalue:
                temp_results[aggregate_levels] = 1
                binning_reads[query_name][current_read_name].append(
                    [temp_mi, temp_ml, current_evalue, aggregate_levels]
                )

            previous_read_name = current_read_name

        update_results(results, sample_index, temp_results, normalise, number_samples)

    if delete_alignments and Path(alignment).exists():
        os.remove(alignment)

    return results, binning_reads
