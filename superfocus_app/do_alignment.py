# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import subprocess
import sys
import time

from pathlib import Path
from collections import defaultdict

# increase CSV load limit
csv.field_size_limit(100000000)


def normalise_counts(data):
    """Normalise query counts based on number of hits.

    Args:
        data (collections.defaultdict): Data to be normalised

    Returns:
        dict: Normalised data

    """
    sum_hits = sum(data.values())

    return {level: data[level] / sum_hits for level in data}


def update_results(results, sample_index, data, normalise, number_samples):
    """Update final result dict.

    Args:
        results (collections.defaultdict): Results that is updated for each sample.
        sample_index (int): Sample index in the result.
        data (dict): Data to be added into the final results.
        normalise (int): 0/1 Normalise or not the data.
        number_samples (int): Number of samples in the analyses.

    Returns:
        collections.defaultdict: Updated results.

    """
    # if any hit was recorded
    if data:
        if normalise and len(data) > 1:
            data = normalise_counts(data)

        # update final result
        for level in data:
            if level not in results:
                results[level] = [0] * number_samples
            results[level][sample_index] += data[level]

    return results


def align_reads(query, output_dir, aligner, database, evalue, threads, fast_mode, WORK_DIRECTORY, amino_acid,
                temp_folder, latency_delay=0):
    """Align FAST(A/Q) file to database.

    Args:
        query (PosixPath): Path to query in FAST(A/Q) file.
        output_dir (PosixPath): Path to alignment output.
        aligner (str): Aligner name.
        database (str): Database name.
        evalue (str): E-value.
        threads (str): Number of threads (default = all).
        fast_mode (str): Fast or sensitive mode (default = 1).
        WORK_DIRECTORY (str): Path to directory where works happens.
        amino_acid (str): 0 input nucleotides, 1 amino acid.
        temp_folder: a temporary directory to write to
        latency_delay: a time delay we will pause between writing and reading the files. This allows a cluster (or NFS mount) to catch up!

    Returns:
        str: Path to alignment that was written.

    """
    # prepare variables
    output_name = "{}/{}_alignments".format(output_dir, query.parts[-1])
    database = "{}_clusters".format(database)
    blast_mode = 'blastp' if amino_acid == '1' else 'blastx'

    if aligner == "diamond":
        #mode_diamond = "" if fast_mode == "1" else "--sensitive"
        database_diamond = "{}/db/static/diamond/{}.db".format(WORK_DIRECTORY, database)

        # align
        # os.system("diamond {} -t {} -d {} -q {} -a {} -p {} -e {} {}".format(blast_mode, temp_folder, database_diamond,
        #                                                                     query, output_name, threads, evalue,
        #                                                                     mode_diamond))
        diamond_blast = [
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
            diamond_blast.append("--sensitive")
        try:
            retcode = subprocess.call(diamond_blast)
            if retcode != 0:
                print("Diamond blast was terminated by signal", retcode, file=sys.stderr)
                sys.exit(retcode)
        except OSError as e:
            print("Diamond blast execution failed:", e, file=sys.stderr)
            sys.exit()

        # note that this is no longer a two step process
        # the daa format is legacy, and diamond supports native
        # tab separated output
        # dump
        #os.system("diamond view -a {}.daa -o {}.m8 -t {} -p {}".format(output_name, output_name, temp_folder, threads))

        # add aligner extension to output
        output_name = "{}.m8".format(output_name)

    elif aligner == 'mmseqs':
        database_mmseqs = f"{WORK_DIRECTORY}/db/static/mmseqs2/{database}.db"
        output_name = f"{output_name}.m8"
        mmseqs_blast = [
            "mmseqs", "easy-search",
            query,
            database_mmseqs,
            output_name,
            temp_folder,
            "--threads", threads,
            "-e", evalue,
        ]
        if fast_mode == "1":
            mmseqs_blast += ["-s", "1.0"]
        try:
            retcode = subprocess.call(mmseqs_blast)
            if retcode != 0:
                print("mmseqs2 blast was terminated by signal", retcode, file=sys.stderr)
                sys.exit(retcode)
        except OSError as e:
            print("mmseqs2 blast execution failed:", e, file=sys.stderr)
            sys.exit()
    elif aligner == "rapsearch":
        mode_rapsearch = "T" if fast_mode == "1" else "F"
        database_rapsearch = "{}/db/static/rapsearch2/{}.db".format(WORK_DIRECTORY, database)

        os.system('rapsearch -a {} -q {} -d {} -o {} -v 250 -z {} -e {} -b 0 -s f'.format(mode_rapsearch, query,
                                                                                          database_rapsearch,
                                                                                          output_name, threads, evalue))
        # add aligner extension to output
        output_name = "{}.m8".format(output_name)
    elif aligner == "blast":
        database_blast = "{}/db/static/blast/{}.db".format(WORK_DIRECTORY, database)

        os.system('{} -db {} -query {} -out {} -outfmt 6 -evalue {} -max_target_seqs 250 -num_threads {}'.format(
            blast_mode, database_blast, query, output_name, evalue, threads))
    else:
        sys.stderr.write(f"FATAL: Aligner {aligner} not known\n")
        sys.exit(0)

    return output_name


def parse_alignments(alignment, results, normalise, number_samples, sample_index, minimum_identity, minimum_alignment,
                     subsystems_translation, aligner, binning_reads, query_name, delete_alignments):
    """Parses alignment.

    Args:
        alignment (PosixPath): Path to sample alignment.
        results (collections.defaultdict): Results that is updated for each sample.
        normalise (int): 0/1 Normalise or not the data.
        number_samples (int): Number of samples in the analyses.
        sample_index (int): Sample index in the result.
        minimum_identity (int): Minimum identity to consider a hit.
        minimum_alignment (int) Minimum alignment (bp) to be consider a hit.
        subsystems_translation (dict): subsystems translation lookup table.
        aligner (str): aligner name.
        binning_reads (collections.defaultdict): reads binning results.
        query_name (str): fasta/q file name used in the alignment.
        delete_alignments (bool): True if files of alignments should be deleted.

    Returns:
        collections.defaultdict: Updated results.

    """
    previous_read_name = None
    best_evalue = None
    temp_results = defaultdict()

    # test for an empty file
    if not os.path.exists(alignment) or os.stat(alignment).st_size == 0:
        return defaultdict(int)

    with open(alignment, encoding='ISO-8859-1') as alignment_file:
        alignment_reader = csv.reader(alignment_file, delimiter='\t')
        if aligner == "rapsearch":
            # skip header
            [next(alignment_reader, None) for _ in range(5)]

        temp_results = defaultdict(int)
        for row in alignment_reader:
            # extract need info from hit
            current_read_name = row[0]
            temp_mi = float(row[2])
            if aligner == 'mmseqs' and temp_mi <= 1 and temp_mi >=0:
                # mmseqs2 reports fraction identity not percent identity!
                temp_mi *= 100
            temp_ml = float(row[3])
            current_evalue = row[10]

            # create needed info
            current_hit = row[1].split("__")
            current_subsystem_id = current_hit[1]  # Subsystem PK on SF database
            current_function_name = current_hit[-1].replace("\n", "").replace("\r", "")
            aggregate_levels = subsystems_translation[current_subsystem_id] + "\t" + current_function_name

            # found a different read
            if current_read_name != previous_read_name:
                if previous_read_name:
                    update_results(results, sample_index, temp_results, normalise, number_samples)
                temp_results = defaultdict(int)
                best_evalue = current_evalue

            # check if alignment wanted
            if temp_mi >= minimum_identity and temp_ml >= minimum_alignment and current_evalue == best_evalue:
                temp_results[aggregate_levels] = 1
                binning_reads[query_name][current_read_name].append([temp_mi, temp_ml, current_evalue,
                                                                     aggregate_levels])

            previous_read_name = current_read_name

        # last group of reads
        update_results(results, sample_index, temp_results, normalise, number_samples)

    # delete file if wanted and it exists
    if delete_alignments and Path(alignment).exists():
        os.remove(alignment)

    return results, binning_reads
