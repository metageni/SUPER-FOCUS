#!/usr/bin/env python3
"""Pipeline orchestration for SUPER-FOCUS functional profiling."""

import logging
import os
import sys
import shutil
import tempfile

from pathlib import Path
from shutil import which
from collections import defaultdict

from superfocus_app.do_alignment import (align_reads,
                                         parse_alignments)
from superfocus_app.utils import (get_subsystems,
                                  is_valid_number,
                                  is_wanted_file,
                                  write_binning,
                                  write_results,
                                  subsample_reads,
                                  aggregate_level,
                                  get_denominators,
                                  add_relative_abundance)

logger = logging.getLogger(__name__)


def _resolve_work_directory(alternate_directory):
    """Resolve the SUPER-FOCUS work directory from flag, env var, or package default.

    Args:
        alternate_directory (str): Value of the --alternate_directory flag (may be empty).

    Returns:
        Path: Resolved work directory.
    """
    if alternate_directory:
        return Path(alternate_directory)
    if 'SUPERFOCUS_DB' in os.environ:
        return Path(os.environ['SUPERFOCUS_DB'])
    return Path(__file__).parent


def _collect_query_files(query_paths):
    """Expand query paths (files and/or directories) into a sorted list of sequence files.

    Args:
        query_paths (list): List of file/directory path strings.

    Returns:
        list: Sorted list of Path objects for valid sequence files.
    """
    files = []
    for f in query_paths:
        p = Path(f)
        if p.is_dir():
            files += [Path(p, x) for x in os.listdir(p)]
        elif p.is_file():
            files.append(p)
    return is_wanted_file(files)


def _validate(aligner, amino_acid, database, threads, evalue, work_directory, run_focus):
    """Validate pipeline parameters before execution.

    Args:
        aligner (str): Aligner name.
        amino_acid (str): '0' or '1'.
        database (str): Cluster size string (e.g. '90').
        threads (str): Thread count or 'all'.
        evalue (str): E-value string.
        work_directory (Path): Resolved work directory.
        run_focus (str): '0' or '1'.

    Returns:
        str or None: Error message if invalid, else None.
    """
    if run_focus != '0':
        return "FOCUS: not available in this version. See https://github.com/metageni/FOCUS"
    if database not in ("90", "95", "98", "100"):
        return f"DATABASE: DB_{database} not valid. Choose DB_90/95/98/100"
    if aligner not in ("diamond", "mmseqs"):
        return f"ALIGNER: {aligner} is not valid. Choose diamond or mmseqs2"
    if not which(aligner):
        return f"ALIGNER: {aligner} is not in PATH"
    if not work_directory.exists():
        return f"WORK_DIRECTORY: {work_directory} does not exist"
    if threads != "all" and not is_valid_number(threads):
        return f"THREADS: {threads} is not a valid number"
    if not is_valid_number(evalue):
        return f"E-VALUE: {evalue} is not valid"
    return None


def run(query_paths, output_directory, prefix, aligner, database, evalue, threads,
        fast_mode, amino_acid, minimum_identity, minimum_alignment, normalise_output,
        alternate_directory, delete_alignments, subsample, latency_wait, run_focus,
        temp_directory):
    """Run the full SUPER-FOCUS functional profiling pipeline.

    Args:
        query_paths (list): List of query file/directory path strings.
        output_directory (str): Path to output directory.
        prefix (str): Output file prefix.
        aligner (str): Aligner name (diamond or mmseqs2).
        database (str): Database cluster size (DB_90, DB_95, DB_98, DB_100).
        evalue (str): E-value threshold.
        threads (str): Thread count or 'all'.
        fast_mode (str): '1' for fast, '0' for sensitive.
        amino_acid (str): '1' for amino acid input, '0' for nucleotide.
        minimum_identity (float): Minimum percent identity for a hit.
        minimum_alignment (int): Minimum alignment length for a hit.
        normalise_output (int): 1 to normalise per-read counts, 0 for raw.
        alternate_directory (str): Alternate database directory path.
        delete_alignments (bool): Delete alignment files after parsing.
        subsample (int or None): Subsample to this many reads, or None.
        latency_wait (int): Seconds to wait after alignment (NFS latency).
        run_focus (str): '0' to skip FOCUS reduction.
        temp_directory (str or None): Explicit temp directory, or None.

    Returns:
        None
    """
    if aligner == 'mmseqs2':
        aligner = 'mmseqs'

    db_size = database.split("_")[-1]
    work_directory = _resolve_work_directory(alternate_directory)
    output_directory = Path(output_directory)

    error = _validate(aligner, amino_acid, db_size, threads, evalue, work_directory, run_focus)
    if error:
        logger.critical(error)
        sys.exit(1)

    output_directory.mkdir(parents=True, exist_ok=True)

    query_files = _collect_query_files(query_paths)
    if not query_files:
        logger.critical("QUERY: no fasta/fna/fastq files found in %s", query_paths)
        sys.exit(1)

    tmp_base = temp_directory or os.environ.get('TMPDIR', '/tmp')
    os.makedirs(tmp_base, exist_ok=True)
    tmpdir = tempfile.mkdtemp(dir=tmp_base)
    logger.info("Using %s as temporary directory", tmpdir)

    try:
        results = defaultdict(list)
        binning_reads = defaultdict(lambda: defaultdict(list))
        subsystems_translation = get_subsystems(Path(work_directory, "db/database_PKs.txt"))

        for counter, temp_query in enumerate(query_files):
            logger.info("1.%d) Working on: %s", counter + 1, temp_query)
            original_query = temp_query
            if subsample:
                temp_query = subsample_reads(temp_query, tmpdir, subsample)

            alignment_name = align_reads(
                temp_query,
                output_dir=output_directory,
                aligner=aligner,
                database=db_size,
                evalue=evalue,
                threads=threads,
                fast_mode=fast_mode,
                WORK_DIRECTORY=work_directory,
                amino_acid=amino_acid,
                temp_folder=tmpdir,
                latency_delay=latency_wait,
            )
            logger.info("   Parsing alignments")
            sample_position = query_files.index(original_query)
            results, binning_reads = parse_alignments(
                alignment_name, results, normalise_output, len(query_files),
                sample_position, minimum_identity, minimum_alignment,
                subsystems_translation, aligner, binning_reads,
                original_query, delete_alignments,
            )

        normalizer = get_denominators(results)
        header_files = query_files + ["{} %".format(x) for x in query_files]
        logger.info("Writing results to %s", output_directory)

        write_binning(binning_reads,
                      "{}/{}binning.xls".format(output_directory, prefix),
                      query_paths, db_size, aligner)

        for level in (1, 2, 3):
            temp_results = aggregate_level(results, level - 1, normalizer)
            write_results(temp_results,
                          ["Subsystem {}".format(level)] + header_files,
                          "{}/{}subsystem_level_{}.xls".format(output_directory, prefix, level),
                          query_paths, db_size, aligner)

        temp_results = add_relative_abundance(results, normalizer)
        write_results(temp_results,
                      ["Subsystem Level 1", "Subsystem Level 2", "Subsystem Level 3", "Function"] + header_files,
                      "{}/{}all_levels_and_function.xls".format(output_directory, prefix),
                      query_paths, db_size, aligner)

    finally:
        shutil.rmtree(tmpdir)

    logger.info("Done")
